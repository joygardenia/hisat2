/*
 * Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
 *
 * This file is part of HISAT-3N.
 *
 * HISAT-3N is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT-3N is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include "alignment_3n_table.h"
#include "utility_3n_table.h"
#include <cassert>
#include <condition_variable>
#include <cstdlib> // size_t
#include <cstring> // memchr, memcpy, memset 等
#include <fcntl.h>
#include <fstream>
#include <iostream> // <-- 用于 cout / cerr
#include <mutex>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;
extern int skipSize;

/**
 * store unique information for one base information with readID, and the
 * quality.
 */
class uniqueID {
public:
  unsigned long long readNameID;
  bool isConverted;
  char quality;
  bool removed;

  uniqueID(unsigned long long InReadNameID, bool InIsConverted, char &InQual) {
    readNameID = InReadNameID;
    isConverted = InIsConverted;
    quality = InQual;
    removed = false;
  }
};

/**
 * basic class to store reference position information
 */
class Position {
  mutex mutex_;

public:
  string chromosome;         // reference chromosome name
  long long int location;    // 1-based position
  char strand;               // +(REF) or -(REF-RC)  正负链
  string convertedQualities; // each char is a mapping quality on this position
                             // for converted base.
  string unconvertedQualities; // each char is a mapping quality on this
                               // position for unconverted base.
  vector<uniqueID>
      uniqueIDs; // each value represent a readName which contributed the base
                 // information. readNameIDs is to make sure no read contribute
                 // 2 times in same position.

  void initialize() {
    chromosome.clear();
    location = -1;
    strand = '?';
    convertedQualities.clear();
    unconvertedQualities.clear();
    vector<uniqueID>().swap(uniqueIDs);
  }

  Position() { initialize(); };

  /**
   * return true if there is mapping information in this reference position.
   */
  bool empty() {
    return convertedQualities.empty() && unconvertedQualities.empty();
  }

  /**
   * set the chromosome, location (position), and strand information.
   */

  void set(string &inputChr, long long int inputLoc) {
    chromosome = inputChr;
    location = inputLoc + 1;
  }

  void set(char inputStrand) { strand = inputStrand; }

  /**
   * binary search of readNameID in readNameIDs.
   * always return a index.
   * if cannot find, return the index which has bigger value than input
   * readNameID.
   */
  int searchReadNameID(unsigned long long &readNameID, int start, int end) {
    if (uniqueIDs.empty()) {
      return 0;
    }
    if (start <= end) {
      int middle = (start + end) / 2;
      if (uniqueIDs[middle].readNameID == readNameID) {
        return middle;
      }
      if (uniqueIDs[middle].readNameID > readNameID) {
        return searchReadNameID(readNameID, start, middle - 1);
      }
      return searchReadNameID(readNameID, middle + 1, end);
    }
    return start; // return the bigger one
  }

  /**
   * with a input readNameID, add it into readNameIDs.
   * if the input readNameID already exist in readNameIDs, return false.
   */
  bool appendReadNameID(PosQuality &InBase, Alignment &InAlignment) {
    int idCount = uniqueIDs.size();
    if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID) {
      uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted,
                             InBase.qual);
      return true;
    }
    int index = searchReadNameID(InAlignment.readNameID, 0, idCount);
    if (uniqueIDs[index].readNameID == InAlignment.readNameID) {
      // if the new base is consistent with exist base's conversion status,
      // ignore otherwise, delete the exist conversion status
      if (uniqueIDs[index].removed) {
        return false;
      }
      if (uniqueIDs[index].isConverted != InBase.converted) {
        uniqueIDs[index].removed = true;
        if (uniqueIDs[index].isConverted) {
          for (int i = 0; i < convertedQualities.size(); i++) {
            if (convertedQualities[i] == InBase.qual) {
              convertedQualities.erase(convertedQualities.begin() + i);
              return false;
            }
          }
        } else {
          for (int i = 0; i < unconvertedQualities.size(); i++) {
            if (unconvertedQualities[i] == InBase.qual) {
              unconvertedQualities.erase(unconvertedQualities.begin() + i);
              return false;
            }
          }
        }
      }
      return false;
    } else {
      uniqueIDs.emplace(uniqueIDs.begin() + index, InAlignment.readNameID,
                        InBase.converted, InBase.qual);
      return true;
    }
  }

  /**
   * append the SAM information into this position.
   */
  void appendBase(PosQuality &input, Alignment &a) {
    mutex_.lock();
    if (appendReadNameID(input, a)) {
      if (input.converted) {
        convertedQualities += input.qual;
      } else {
        unconvertedQualities += input.qual;
      }
    }
    mutex_.unlock();
  }
};

/**
 * store all reference position in this class.
 */
class Positions {
public:
  vector<Position *>
      refPositions;       // the pool of all current reference position.
  string chromosome;      // current reference chromosome name.
  string cout_Chromosome; // cout chromosome name
  long long int
      location; // current location (position) in reference chromosome.
  char lastBase =
      'X'; // the last base of reference line. this is for CG_only mode.
  SafeQueue<string *> linePool;     // pool to store unprocessed SAM line.
  SafeQueue<string *> freeLinePool; // pool to store free string pointer for SAM
                                    // line.    一行 mode:line
  SafeQueue<Alignment *> freeAlignmentPool; // new
  streampos startPos;                       // new
  vector<Position> group1;
  vector<Position> group2;
  Position *ref;
  LRUCache LRU;
  SafeQueue<Position *>
      freePositionPool; // pool to store free position pointer mode:refPosition
                        // 一个位置 for reference position.
  SafeQueue<Position *>
      outputPositionPool; // pool to store the reference position which is
                          // 是防止多线程删除出事吗 loaded and ready to output.

  bool working;
  mutex mutex_;
  long long int refCoveredPosition; // this is the last position in reference
                                    // chromosome we loaded in refPositions.
  ifstream refFile;
  condition_variable refcv;
  condition_variable linecv;
  vector<mutex *> workerLock; // one lock for one worker thread.
  int nThreads = 1;
  ChromosomeFilePositions
      chromosomePos; // store the chromosome name and it's streamPos. To quickly
                     // find new chromosome in file.
  bool addedChrName = false;
  bool removedChrName = false;

  Positions(int inputNThreads, bool inputAddedChrName, bool inputRemovedChrName)
      : group1(6000), group2(6000), LRU(16) {
    ref = group1.data();
    working = true;
    nThreads = inputNThreads;
    addedChrName = inputAddedChrName;
    removedChrName = inputRemovedChrName;
  }
  ~Positions() {
    string *line;
    while (freeLinePool.popFront(line)) {
      delete line;
    }
    Alignment *newAlignment;
    while (freeAlignmentPool.popFront(newAlignment)) {
      delete newAlignment;
    }
  }
  /**
   * given the target Position output the corresponding position index in
   * refPositions.给定一个目标坐标（targetPos），计算它相对于当前参考片段起点的偏移量（index）
   */
  int getIndex(long long int &targetPos) {
    int firstPos = ref->location;
    return targetPos - firstPos;
  }

  void LoadChromosomeNamesPos(const string &refFileName) {
    if (!readFaiIfFresh(refFileName)) {
      scanAndBuildFai(refFileName);
    }
    for (auto &c : chromosomePos.pos) {
      if (loadingBlockSize % c.lineBases != 0) {
        cerr << "warning: page size " << loadingBlockSize
             << " is not a multiple of lineBases " << c.lineBases
             << " (chromosome " << c.chromosome << "); pages may misalign."
             << endl;
        break;
      }
    }
  }

  /**
   * 读取 samtools 风格的 <ref>.fai。存在且不早于参考基因组时返回 true。
   * 格式:name \t length \t offset \t linebases \t linewidth
   */
  bool readFaiIfFresh(const string &refFileName) {
    string faiName = refFileName + ".fai";
    struct stat refStat, faiStat;
    if (stat(refFileName.c_str(), &refStat) != 0)
      return false;
    if (stat(faiName.c_str(), &faiStat) != 0)
      return false;
    if (faiStat.st_mtime < refStat.st_mtime)
      return false;

    ifstream fai(faiName);
    if (!fai.good())
      return false;

    string line;
    while (getline(fai, line)) {
      if (line.empty() || line.front() == '#')
        continue;
      size_t t1 = line.find('\t');
      if (t1 == string::npos)
        continue;
      size_t t2 = line.find('\t', t1 + 1);
      if (t2 == string::npos)
        continue;
      size_t t3 = line.find('\t', t2 + 1);
      if (t3 == string::npos)
        continue;
      size_t t4 = line.find('\t', t3 + 1);
      if (t4 == string::npos)
        continue;

      string chr = line.substr(0, t1);
      streampos off = (streampos)stoll(line.substr(t2 + 1, t3 - t2 - 1));
      int lb = stoi(line.substr(t3 + 1, t4 - t3 - 1));
      int lw = stoi(line.substr(t4 + 1));
      if (lb <= 0 || lw <= 0)
        continue;
      chromosomePos.append(chr, off, lb, lw);
    }
    if (chromosomePos.pos.empty())
      return false;
    chromosomePos.sort();
    chromosome.clear();
    return true;
  }

  void scanAndBuildFai(const string &refFileName) {
    int fd = open(refFileName.c_str(), O_RDONLY);
    if (fd < 0)
      throw runtime_error("open failed");

    struct stat sb;
    fstat(fd, &sb);
    size_t size = sb.st_size;

    char *data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (data == MAP_FAILED)
      throw runtime_error("mmap failed");

    char *end = data + size;
    char *p = data;

    const size_t WINDOW = 8 * 1024 * 1024;
    const size_t PAGE = (size_t)sysconf(_SC_PAGESIZE);
    char *releasedUpTo = data;
    auto releaseBefore = [&](char *upto) {
      char *b = (char *)((size_t)upto & ~(PAGE - 1));
      if (b > releasedUpTo) {
        madvise(releasedUpTo, b - releasedUpTo, MADV_DONTNEED);
        releasedUpTo = b;
      }
    };
    auto findGt = [&](char *from) -> char * {
      char *q = from;
      while (q < end) {
        size_t win = (size_t)(end - q) < WINDOW ? (size_t)(end - q) : WINDOW;
        char *g = (char *)memchr(q, '>', win);
        if (g)
          return g;
        q += win;
        releaseBefore(q);
      }
      return NULL;
    };

    ofstream fai(refFileName + ".fai");

    while (p < end) {
      char *gt = findGt(p);
      if (gt == NULL)
        break;

      if (gt != data && gt[-1] != '\n') {
        p = gt + 1;
        continue;
      }
      char *nl = (char *)memchr(gt, '\n', end - gt);
      if (!nl)
        nl = end;

      string header(gt + 1, nl);
      string chr;
      for (char c : header) {
        if (isspace((unsigned char)c))
          break;
        chr += c;
      }

      char *lineStart = nl + 1;
      char *lineEnd = (char *)memchr(lineStart, '\n', end - lineStart);
      while (lineEnd != NULL && lineEnd == lineStart) {
        lineStart = lineEnd + 1;
        lineEnd = (char *)memchr(lineStart, '\n', end - lineStart);
      }
      if (lineEnd == NULL)
        lineEnd = end;

      int raw = (int)(lineEnd - lineStart);
      int cr = (raw > 0 && lineStart[raw - 1] == '\r') ? 1 : 0;
      int lb = raw - cr;
      int lw = raw + 1; // 相邻两行起点的字节距离,自动兼容 LF / CRLF

      char *nextGt = findGt(lineStart);
      char *limit = nextGt ? nextGt : end;

      if (lb > 0) {
        streampos off = lineStart - data;
        // len = 整行数 * lb + 末行碱基数,由字节数 O(1) 推出,无需再扫一遍
        long long bytes = limit - lineStart;
        long long rem = bytes % lw;
        long long tail = rem > 0 ? rem - 1 - cr : 0;
        if (tail < 0)
          tail = 0;
        long long len = (bytes / lw) * lb + tail;

        chromosomePos.append(chr, off, lb, lw);
        if (fai.good())
          fai << chr << '\t' << len << '\t' << (long long)off << '\t' << lb
              << '\t' << lw << '\n';
      }

      p = limit;
      releaseBefore(p);
    }

    munmap(data, size);
    chromosomePos.sort();
    chromosome.clear();
  }
  /**
   * get a fasta line (not header), append the bases to positions.
   */
  void appendRefPosition(string &line, int &index_group) {
    Position *newPos;
    // check the base one by one
    char *b;
    for (int i = 0; i < line.size(); i++) {
      // getFreePosition(newPos);
      // newPos = new Position(); // TODO:修改成freepositionpool.popfront
      newPos = ref + (index_group++);
      newPos->set(cout_Chromosome,
                  location + i); // 在执行loadstart时location是固定的位置
      b = &line[i];
      if (CG_only) {
        if (lastBase == 'C' && *b == 'G') {
          (ref + (index_group - 2))
              ->set('+'); // 修改上一个碱基（刚刚存进去的那个位置）的属性
          newPos->set('-');
        }
      } else {
        if (*b == convertFrom) {
          newPos->set('+'); // 正链：A、T、C、G 的原始方向
        } else if (*b == convertFromComplement) {
          newPos->set('-');
        }
      }
      // refPositions.push_back( // 创建2行120个refposition
      //     newPos); // 把 newPos（一个 Position* 指针）加到 vector 尾部
      lastBase = *b;
    }
    location += line.size(); // ref的块读入位移，不能乱改
  }

  /**
   * if we can go through all the workerLock, that means no worker is
   * appending new position.
   * 有道理但不多，这eee牺牲性能实现正确性，这不该是我做的事情吗
   */
  void appendingFinished() {
    for (int i = 0; i < nThreads; i++) {
      workerLock[i]->lock();
      workerLock[i]->unlock();
    }
  }

  /**
   * the output function for output thread.
   */
  void outputFunction(string outputFileName) {
    ostream *out_ = &cout;
    out_ = &cout;
    ofstream tableFile;
    if (!outputFileName.empty()) {
      tableFile.open(outputFileName, ios_base::out);
      out_ = &tableFile;
    }

    *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunc"
             "onvertedBaseQualities\tunconvertedBaseCount\n";
    Position *pos;
    while (working) {
      if (outputPositionPool.popFront(pos)) {
        *out_ << pos->chromosome << '\t' << to_string(pos->location) << '\t'
              << pos->strand << '\t' << pos->convertedQualities << '\t'
              << to_string(pos->convertedQualities.size()) << '\t'
              << pos->unconvertedQualities << '\t'
              << to_string(pos->unconvertedQualities.size()) << '\n';
        returnPosition(pos);
      } else {
        this_thread::sleep_for(std::chrono::microseconds(1));
      }
    }
    tableFile.close();
  }

  void writeHeader(string outputFileName) {
    if (outputFileName.empty())
      return;
    ofstream tableFile(outputFileName, ios::out); // 覆盖模式，文件会被清空
    tableFile
        << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\t"
           "unconvertedBaseQualities\tunconvertedBaseCount\n";

    tableFile.close();
  }

  void _moveAllToprint(string outputFileName) { //&
    ostream *out_ = &cout;

    ofstream tableFile;
    if (!outputFileName.empty()) {
      tableFile.open(outputFileName, ios::app); // append追加模式
      out_ = &tableFile;
    }
    int blockstart = 0;
    int blockend = loadingBlockSize - 1;
    //    int blockstart =
    //        (newAlignment.location % loadingBlockSize - 1 +
    //        loadingBlockSize) % loadingBlockSize; //
    //        原始是front_letter，后一页不方便得到
    //    int blockend =
    //        blockstart + newAlignment.endRefPos - newAlignment.startRefPos;
    Position *pos;
    //*out_ << "size of refposition : " << refPositions.size() << "\n";
    //*out_ << "blockstart : " << blockstart << "\n";
    //*out_ << "blockend : " << blockend << "\n";
    //*out_ << "newAlignment.letter : " << newAlignment.letter << "\n";
    for (int index = blockstart; index <= blockend; index++) {
      pos = ref + index;
      if (pos->empty() || pos->strand == '?') {
        continue;
      } else {
        // vector<uniqueID>().swap(
        //     group[index].uniqueIDs); // 清空 vector 并释放它占用的堆内存”
        //  的高效写法，比 clear() 更彻底。

        *out_ << pos->chromosome << '\t' << to_string(pos->location) << '\t'
              << pos->strand << '\t' << pos->convertedQualities << '\t'
              << to_string(pos->convertedQualities.size()) << '\t'
              << pos->unconvertedQualities << '\t'
              << to_string(pos->unconvertedQualities.size()) << '\n';
      }
    }
    for (int i = 0; i < loadingBlockSize; i++)
      (ref + i)->initialize();
    tableFile.close();
    // refPositions.clear();
  }
  /**
   * move the position which position smaller than refCoveredPosition -
   * loadingBlockSize, output it.
   */
  void moveBlockToOutput() {
    if (refPositions.empty()) {
      return;
    }
    int index;
    for (index = 0; index < refPositions.size(); index++) {
      if (refPositions[index]->location <
          refCoveredPosition - loadingBlockSize) {
        if (refPositions[index]->empty() ||
            refPositions[index]->strand == '?') {
          returnPosition(refPositions[index]);
        } else {
          outputPositionPool.push(refPositions[index]);
        }
      } else {
        break;
      }
    }
    if (index != 0) {
      refPositions.erase(refPositions.begin(), refPositions.begin() + index);
    }
  }

  /**
   * move all the refPosition into output pool.
   */
  void moveAllToOutput() {
    if (refPositions.empty()) {
      return;
    }
    for (int index = 0; index < refPositions.size(); index++) {
      if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
        returnPosition(refPositions[index]);
      } else {
        vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
        outputPositionPool.push(refPositions[index]);
      }
    }
    refPositions.clear();
  }

  void loadTestChromosome(string targetChromosome) {
    refFile.clear();
    startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
    cout_Chromosome = targetChromosome;
    // chromosome = targetChromosome; // 可以去掉，暂时别管
  }
  void startload(long long int spage) {
    long long int pageStart = loadingBlockSize * spage;
    refFile.seekg(chromosomePos.byteOffsetOf(cout_Chromosome, pageStart),
                  ios::beg);
    string line;
    lastBase = 'X';
    location = pageStart;

    int index_group = 0;
    int CountLine = loadingBlockSize / chromosomePos.lineBasesOf(cout_Chromosome);
    while (refFile.good() && CountLine) {
      getline(refFile, line);
      if (!line.empty() && line.back() == '\r')
        line.pop_back();
      if (line.empty()) {
        continue;
      }
      if (line.front() == '>') { // meet next chromosome, return.
        return;
      }
      // change all base to upper case
      for (int i = 0; i < line.size(); i++) {
        line[i] = toupper(line[i]);
      }
      appendRefPosition(line, index_group);
      CountLine--;
    }
  }

  /**
   * 把该读段落在本页的碱基切片累加进参考位点。
   * 绝对坐标 X = location + refPos (1-based) -> 本页槽位 slot = X - 1 - page*B。
   */
  void appendPositions(const PageSlice &slice, int page) {
    Alignment &a = *slice.a;
    if (!a.mapped || a.bases.empty()) {
      return;
    }
    long long pageStart0 = (long long)page * loadingBlockSize;
    for (int i = slice.lo; i < slice.hi; i++) {
      PosQuality *b = &a.bases[i];
      if (b->remove) {
        continue;
      }
      int slot = (int)(a.location + b->refPos - 1 - pageStart0);
      assert(slot >= 0 && slot < loadingBlockSize);
      Position *pos = ref + slot;

      if (pos->strand == '?') {
        // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
        continue;
      }
      pos->appendBase(*b, a); // 把这条 read
                              // 对这个位置的观测值（碱基种类、测序质量、甲基化标记等）累积到
                              // Position 对象里。后续统计如 coverage、甲基化比例等都是在这些
                              // 累计信息上进行。
    }
  }

  /**
   * get a string pointer from freeLinePool, if freeLinePool is empty, make a
   * new string pointer.
   */
  void getFreeStringPointer(string *&newLine) {
    if (freeLinePool.popFront(newLine)) {
      return;
    } else {
      newLine = new string();
    }
  }

  /**
   * get a Position pointer from freePositionPool, if freePositionPool is
   * empty, make a new Position pointer.
   */
  void getFreePosition(Position *&newPosition) {
    while (outputPositionPool.size() >= 10000) {
      this_thread::sleep_for(std::chrono::microseconds(1));
    }
    if (freePositionPool.popFront(newPosition)) {
      return;
    } else {
      newPosition = new Position();
    }
  }

  /**
   * return the line to freeLinePool
   */
  void returnLine(string *line) {
    line->clear();
    freeLinePool.push(line);
  }

  /**
   * return the position to freePositionPool.
   */
  void returnPosition(Position *pos) {
    pos->initialize();
    freePositionPool.push(pos);
  }

  void
  getFreeAlignment(Alignment *&newAlignment) { // 对一个指向 Alignment
                                               // 的指针的引用,只用*无法改变外部
    if (freeAlignmentPool.popFront(newAlignment)) {
      return;
    } else {
      newAlignment = new Alignment;
    }
  }
  void returnAlignment(Alignment *newAlignment) {
    newAlignment->initialize();
    //: 写alignment的clear函数,已经写好了
    freeAlignmentPool.push(newAlignment);
  }
};

#endif // POSITION_3N_TABLE_H
