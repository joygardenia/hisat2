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

  /**
   * given reference line (start with '>'), extract the chromosome information.
   * this is important when there is space in chromosome name. the SAM
   * information only contain the first word.
   */
  string getChrName(string &inputLine) {
    string name;
    for (int i = 1; i < inputLine.size(); i++) {
      char c = inputLine[i];
      if (isspace(c)) {
        break;
      }
      name += c;
    }

    if (removedChrName) {
      if (name.find("chr") == 0) {
        name = name.substr(3);
      }
    } else if (addedChrName) {
      if (name.find("chr") != 0) {
        name = string("chr") + name;
      }
    }
    return name;
  }

  void _LoadChromosomeNamesPos() {
    string line;
    // while (refFile.good()) {
    while (getline(refFile, line)) {
      // getline(refFile, line);
      // if (line.front() == '>') {
      if (!line.empty() &&
          line.front() == '>') { // this line is chromosome name
        chromosome = getChrName(line);
        streampos currentPos = refFile.tellg();
        chromosomePos.append( // vector
            chromosome,
            currentPos); // File position for char
                         // streams.每条染色体的名字和字节流位置？dump文件中的位置
      }
    }
    chromosomePos.sort();
    chromosome.clear();
  }
  void LoadChromosomeNamesPos(const string &refFileName) {
    int fd = open(refFileName.c_str(), O_RDONLY);
    if (fd < 0)
      throw runtime_error("open failed");

    struct stat sb;
    fstat(fd, &sb);
    size_t size = sb.st_size;

    char *data = (char *)mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);

    char *p = data;
    char *end = data + size;

    /*while (p < end) {

     if (*p == '>') {
        char *header_start = p + 1;

        // 找到换行
        char *header_end =
            (char *)memchr(header_start, '\n', end - header_start);
        if (!header_end)
          header_end = end; // 最后一行无换行也成立

        // 构造头部字符串： header_start → header_end-1
        string header(header_start, header_end);

        // 提取真正的染色体名（停在第一个空白）
        string chr;
        for (char c : header) {
          if (isspace((unsigned char)c))
            break;
          chr += c;
        }

        // 得到序列起始位置（下一行）
        streampos seq_pos = (header_end - data) + 1;

        // 加入表
        chromosomePos.append(chr, seq_pos);
      }

      // 移动到下一行
      char *next = (char *)memchr(p, '\n', end - p);
      if (!next)
        break;
      p = next + 1;
    }*/
    while (p < end) {
      char *next = (char *)memchr(p, '\n', end - p);
      if (!next)
        break;

      if (*p == '>') {
        char *header_start = p + 1;
        string header(header_start, next); // 直接用 next，不用再找一次

        string chr;
        for (char c : header) {
          if (isspace((unsigned char)c))
            break;
          chr += c;
        }
        streampos seq_pos = (next - data) + 1;
        chromosomePos.append(chr, seq_pos);
      }

      p = next + 1;
    }

    chromosomePos.sort();
    chromosome.clear();

    munmap(data, size);
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

  // my masterpiece
  void deleteRefPos() {
    for (auto p : refPositions) {
      delete p; // 释放每个 Position*
    }
    refPositions.clear(); // 清空容器
  }
  void moveAllToprint(string outputFileName, Alignment &newAlignment) {
    ostream *out_ = &cout;

    ofstream tableFile;
    if (!outputFileName.empty()) {
      tableFile.open(outputFileName, ios_base::out);
      out_ = &tableFile;
    }
    int blockstart = newAlignment.front_letter - 1;
    int blockend = newAlignment.front_letter + newAlignment.endRefPos -
                   newAlignment.startRefPos - 1;
    *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunc"
             "onvertedBaseQualities\tunconvertedBaseCount\n";
    Position *pos;
    *out_ << "size of refposition : " << refPositions.size() << "\n";
    *out_ << "blockstart : " << blockstart << "\n";
    *out_ << "blockend : " << blockend << "\n";
    *out_ << "newAlignment.front_letter : " << newAlignment.front_letter
          << "\n";
    for (int index = blockstart; index <= blockend; index++) {
      if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
        continue;
      } else {
        vector<uniqueID>().swap(
            refPositions[index]
                ->uniqueIDs); // 清空 vector 并释放它占用的堆内存”
                              // 的高效写法，比 clear() 更彻底。
        pos = refPositions[index];
        *out_ << pos->chromosome << '\t' << to_string(pos->location) << '\t'
              << pos->strand << '\t' << pos->convertedQualities << '\t'
              << to_string(pos->convertedQualities.size()) << '\t'
              << pos->unconvertedQualities << '\t'
              << to_string(pos->unconvertedQualities.size()) << '\n';
      }
    }
    tableFile.close();
    refPositions.clear();
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
    // refCoveredPosition = 2 * loadingBlockSize;
    refFile.seekg(startPos, ios::beg); // 跳转到该染色体位置
    // refCoveredPosition = 2 * loadingBlockSize;
    string line;
    lastBase = 'X';
    location = loadingBlockSize * spage;

    int index_group = 0;
    int CountLine = loadingBlockSize / 60; // 一页里面多少行
    refFile.seekg(spage * skipSize * CountLine, ios::cur);
    while (refFile.good() && CountLine) {
      getline(refFile, line);
      if (line.front() == '>') { // meet next chromosome, return.
        return;
      } else {
        if (line.empty()) {
          continue;
        }

        // change all base to upper case
        for (int i = 0; i < line.size(); i++) {
          line[i] = toupper(line[i]);
        }
        appendRefPosition(line, index_group);
        CountLine--;
      }
    }
  }

  /**
   * load more Position (loadingBlockSize bp) to positions
   * if we meet next chromosome, return false. Else, return ture.
   */
  void loadMore(int CountLine) {
    // refCoveredPosition += loadingBlockSize;
    string line;
    while (refFile.good() && CountLine) {
      getline(refFile, line);
      if (line.front() == '>') { // meet next chromosome, return.
        return;
      } else {
        if (line.empty()) {
          continue;
        }

        // change all base to upper case
        for (int i = 0; i < line.size(); i++) {
          line[i] = toupper(line[i]);
        }
        // appendRefPosition(line);
        CountLine--;
        //        if (location >=
        //            refCoveredPosition) { //
        //            超过了已读取块的总量，return,重新readmore
        //          return;
        //        }
      }
    }
  }

  /**
   * add position information from Alignment into ref position.
   */
  void appendPositions(Alignment &newAlignment) {
    if (!newAlignment.mapped || newAlignment.bases.empty()) {
      return;
    }

    int start = 0;
    int end = newAlignment.sequence.size();
    if (newAlignment.middleRefPos != -1) { // 存在分页
      // FIXME:判断是前面的页面还是后面的页面，修改start end
      if (newAlignment.isFirstPage) {
        end = newAlignment.middleRefPos;
        newAlignment.endRefPos = newAlignment.bases[end].refPos;
      } else {
        // startPos =
        // (newAlignment.bases[newAlignment.middleRefPos].refPos+newAlignment.location)%loadingBlockSize;//在alignment
        // 的base[i]全部改成以middleRefPos的初始地址为偏移量了
        newAlignment.location =
            newAlignment.bases[newAlignment.middleRefPos].refPos +
            newAlignment.location;
        newAlignment.bases[newAlignment.middleRefPos].refPos = 0;

        start = newAlignment.middleRefPos;
        newAlignment.startRefPos = newAlignment.bases[start].refPos;
      }
    }
    int index = getIndex(newAlignment.location);
    long long int startPos = newAlignment.location; // 1-based position
    // find the first reference position in pool.
    for (int i = start; i < end;
         i++) { //<end 把原来的前面页面的middleRefPos给磨没了
      PosQuality *b = &newAlignment.bases[i];
      if (b->remove) {
        continue;
      }

      Position *pos = ref + (index + b->refPos);
      assert(pos->location == startPos + b->refPos);

      if (pos->strand == '?') {
        // this is for CG-only mode. read has a 'C' or 'G' but not 'CG'.
        continue;
      }
      pos->appendBase(
          newAlignment.bases[i],
          newAlignment); // 把这条 read
                         // 对这个位置的观测值（碱基种类、测序质量、甲基化标记等）累积到
                         // Position 对象里。后续统计如
                         // coverage、甲基化比例等都是在这些累计信息上进行。
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
  /**
   * this is the working function.
   * it take the SAM line from linePool, parse it.
   */
  void append(int threadID) {
    string line;
    Alignment newAlignment;

    //    while (working) {
    //      {
    //        std::unique_lock<std::mutex> lock(*workerLock[threadID]);
    //        // No matching constructor for initialization of
    //        //
    //        'std::unique_lock<std::mutex>'------->>>>>>>>>>>没有找到能匹配的构造函数,即构造函数要求传入的是
    //        // mutex对象的引用，而非一个指针（std::mutex**）
    //        linecv.wait(lock, [&] { return linePool.popFront(line); });
    //      }
    //      {
    //        std::unique_lock<std::mutex> lock(*workerLock[threadID]);
    //        refcv.wait(lock, [&] { return !refPositions.empty(); });
    //        newAlignment.parse(line);
    //        returnLine(line);
    //        appendPositions(newAlignment);
    //      }
    //    }

    newAlignment.parse(&line);
    appendPositions(newAlignment);
  }
};

#endif // POSITION_3N_TABLE_H
