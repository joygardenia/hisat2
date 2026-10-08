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
#include "alignment_3n_table.h"
#include "position_3n_table.h"
#include "utility_3n_table.h"
#include <cxxabi.h>
#include <execinfo.h>
#include <getopt.h>
#include <iostream>
#include <map>
#include <thread>
#include <vector>
using namespace std;

string alignmentFileName;
bool standardInMode = false;
string refFileName;
string outputFileName;
bool uniqueOnly = false;
bool multipleOnly = false;
bool CG_only = false;
char convertFrom = '0';
char convertTo = '0';
char convertFromComplement;
char convertToComplement;
bool addedChrName = false;
bool removedChrName = false;
Positions *positions;
int nThreads = 1;

// 添加代码
long long int loadingBlockSize = 6000;

bool fileExist(string &filename) {
  ifstream file(filename);
  return file.good();
}

enum { ARG_ADDED_CHRNAME = 256, ARG_REMOVED_CHRNAME };

static const char *short_options = "s:r:t:b:umcp:h";
static struct option long_options[]{
    {"alignments", required_argument, 0, 'a'},
    {"ref", required_argument, 0, 'r'},
    {"output-name", required_argument, 0, 'o'},
    {"base-change", required_argument, 0, 'b'},
    {"unique-only", no_argument, 0, 'u'},
    {"multiple-only", no_argument, 0, 'm'},
    {"CG-only", no_argument, 0, 'c'},
    {"threads", required_argument, 0, 'p'},
    {"added-chrname", no_argument, 0, ARG_ADDED_CHRNAME},
    {"removed-chrname", no_argument, 0, ARG_REMOVED_CHRNAME},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}};

static void printHelp(ostream &out) {
  out << "hisat-3n-table developed by Yun (Leo) Zhang" << endl;
  out << "Usage:" << endl
      << "hisat-3n-table [options]* --alignments <alignmentFile> --ref "
         "<refFile> --output-name <outputFile> --base-change <char1,char2>"
      << endl
      << "  <alignmentFile>           SORTED SAM filename. Please enter '-' "
         "for standard input."
      << endl
      << "  <refFile>                 reference file (should be FASTA format)."
      << endl
      << "  <outputFile>              file name to save the 3n table (tsv "
         "format). By default, alignments are written to the “standard out” or "
         "“stdout” filehandle (i.e. the console)."
      << endl
      << "  <chr1,chr2>               the char1 is the nucleotide converted "
         "from, the char2 is the nucleotide converted to."
      << endl;
  out << "Options (defaults in parentheses):" << endl
      << " Input:" << endl
      << "  -u/--unique-only          only count the base which is in unique "
         "mapped reads."
      << endl
      << "  -m/--multiple-only        only count the base which is in multiple "
         "mapped reads."
      << endl
      << "  -c/--CG-only              only count CG and ignore CH in reference."
      << endl
      << "  --added-chrname           please add this option if you use "
         "--add-chrname during HISAT-3N alignment."
      << endl
      << "  --removed-chrname         please add this option if you use "
         "--remove-chrname during HISAT-3N alignment."
      << endl
      << "  -p/--threads <int>        number of threads to launch (1)." << endl
      << "  -h/--help                 print this usage message." << endl;
}

static void parseOption(int next_option, const char *optarg) {
  switch (next_option) {
  case 'a': {
    alignmentFileName = optarg;
    if (alignmentFileName == "-") {
      standardInMode = true; // 从标准输入中读取dump文件
      break;
    }
    if (!fileExist(alignmentFileName)) {
      cerr << "The alignment file is not exist." << endl;
      throw(1);
    }
    break;
  }
  case 'r': {
    refFileName = optarg;
    if (!fileExist(refFileName)) {
      cerr << "reference (FASTA) file is not exist." << endl;
      throw(1);
    }
    break;
  }
  case 'o':
    outputFileName = optarg;
    break;
  case 'b': {
    string arg = optarg;
    if (arg.size() != 3 || arg[1] != ',') {
      cerr << "Error: expected 2 comma-separated "
           << "arguments to --base-change option (e.g. C,T), got " << arg
           << endl;
      throw 1;
    }
    convertFrom = toupper(arg.front());
    convertTo = toupper(arg.back());
    break;
  }
  case 'u': {
    uniqueOnly = true;
    break;
  }
  case 'm': {
    multipleOnly = true;
    break;
  }
  case 'c': {
    CG_only = true;
    break;
  }
  case 'h': {
    printHelp(cerr);
    throw 0;
  }
  case 'p': {
    cerr << "[debug][p] stoi input = '" << optarg << "'" << endl;
    nThreads = stoi(optarg);

    if (nThreads < 1) {
      nThreads = 1;
    }
    break;
  }
  case ARG_ADDED_CHRNAME: {
    addedChrName = true;
    break;
  }
  case ARG_REMOVED_CHRNAME: {
    removedChrName = true;
    break;
  }
  default:
    printHelp(cerr);
    throw 1;
  }
}

static void parseOptions(int argc, const char **argv) {
  int option_index = 0;
  int next_option;
  while (true) {
    next_option = getopt_long(argc, const_cast<char **>(argv), short_options,
                              long_options, &option_index);
    if (next_option == -1)
      break;
    parseOption(next_option, optarg);
  }

  // check filenames
  if (refFileName.empty() || alignmentFileName.empty()) {
    cerr << "No reference or SAM file specified!" << endl;
    printHelp(cerr);
    throw 1;
  }

  // give a warning for CG-only
  if (CG_only) {
    if (convertFrom != 'C' || convertTo != 'T') {
      cerr << "Warning! You are using CG-only mode. The the --base-change "
              "option is set to: C,T"
           << endl;
      convertFrom = 'C';
      convertTo = 'T';
    }
  }

  // check if --base-change is empty
  if (convertFrom == '0' || convertTo == '0') {
    cerr << "the --base-change argument is required." << endl;
    throw 1;
  }

  if (removedChrName && addedChrName) {
    cerr << "Error: --removed-chrname and --added-chrname cannot be used at "
            "the same time"
         << endl;
    throw 1;
  }

  // set complements
  convertFromComplement = asc2dnacomp[convertFrom];
  convertToComplement = asc2dnacomp[convertTo];
}

/**
 * give a SAM line, extract the chromosome and position information.
 * return true if the SAM line is mapped. return false if SAM line is not maped.
 */
bool getSAMChromosomePos(string *line, string &chr, long long int &pos) {
  int startPosition = 0;
  int endPosition = 0;
  int count = 0;

  while ((endPosition = line->find("\t", startPosition)) !=
         string::npos) { // 循环查找每个 \t（制表符），说明要逐列提取字段
    if (count == 2) {
      chr = line->substr(startPosition, endPosition - startPosition);
    } else if (count == 3) {
      pos = stoll(line->substr(startPosition,
                               endPosition - startPosition)); // string->int
      if (chr == "*") {
        return false;
      } else {
        return true;
      }
    }
    startPosition = endPosition + 1;
    count++;
  }
  return false;
}

// 添加代码
string *line;
string samChromosome; // the chromosome name of current SAM line.
long long int
    samPos; // the position of current SAM line.某条染色体某个碱基处开始
void forwardSam(istream *alignmentFile, Positions *positions,
                string *&nextline) {
  int startPosition = 0;
  int endPosition = 0;
  int count = 0;
  positions->getFreeStringPointer(line);
  while (alignmentFile->good()) { // 可以安全读取或写入流
    if (!getline(*alignmentFile, *line)) {
      positions->returnLine(line);
      break;
    }
    if (!line->empty() && line->back() == '\r')
      line->pop_back();
    if (line->empty() || line->front() == '@') {
      continue;
    }
    nextline = line; // 得到第一条read（可识别
    break;
  }
  while ((endPosition = line->find("\t", startPosition)) !=
         string::npos) { // 循环查找每个 \t（制表符），说明要逐列提取字段
    if (count == 2) {
      samChromosome = line->substr(startPosition, endPosition - startPosition);
      break;
    }
    startPosition = endPosition + 1;
    count++;
  }
}
// 添加代码
const int linePoolHighWater = 100000;

void processLinePool(Positions *positions);

bool readNewSamChromosome(istream *alignmentFile, Positions *positions,
                          string *&nextline) {
  // 从 SAM
  // 文件中读取一批有效行，跳过空行和头信息，遇到下一条属于新染色体时就停下，把这一行留给下一轮处理，同时用一个对象池/linePool管理
  // std::string* 的复用
  if (alignmentFile->good() && nextline) {
    positions->linePool.push(nextline);
    nextline = nullptr;
    positions->chromosome = samChromosome;
  }
  while (alignmentFile->good()) { // 可以安全读取或写入流
    if (positions->linePool.size() > linePoolHighWater) {
      // 只排空 linePool。此处不可 popAllNodesFromTail:
      // 最新页随后仍可能被写入,提前刷出会让同一位点分两次输出。
      processLinePool(positions);
    }
    positions->getFreeStringPointer(line);
    if (!getline(*alignmentFile, *line)) {
      positions->returnLine(line);

      return false;
    }
    if (!line->empty() && line->back() == '\r')
      line->pop_back();
    if (line->empty() || line->front() == '@') {
      positions->returnLine(line);
      continue;
    }
    if (!getSAMChromosomePos(line, samChromosome, samPos)) {
      positions->returnLine(line);
      break;
    }
    if (samChromosome != positions->chromosome) { // 结束这条染色体的sam读取
      nextline = line;
      break;
    }
    positions->linePool.push(line);
  }
  return true;
}

// 添加代码
const int outputPoolHighWater = 20000;

void processLinePool(Positions *positions) {
  Alignment *newAlignment;
  while (positions->linePool.popFront(line)) {
    while (positions->LRU.outputPool.size() > outputPoolHighWater) {
      this_thread::sleep_for(std::chrono::microseconds(10));
    }
    positions->getFreeAlignment(newAlignment);
    bool parsed = newAlignment->parse(line); // 解析
    positions->returnLine(line);
    if (!parsed || !newAlignment->mapped || newAlignment->bases.empty()) {
      positions->returnAlignment(newAlignment);
      continue;
    }
    // 按页把读段的碱基切成若干段(不复制读段对象)。kept 碱基的 refPos 单调不减,
    // 故同一页对应连续一段; remove 的碱基落在哪一侧都会被跳过。
    long long loc = newAlignment->location;
    int nb = (int)newAlignment->bases.size();
    int nSlices = 0, cur = -1;
    for (int i = 0; i < nb; ++i) {
      if (newAlignment->bases[i].remove)
        continue;
      int p =
          (int)((loc + newAlignment->bases[i].refPos - 1) / loadingBlockSize);
      if (p != cur) {
        ++nSlices;
        cur = p;
      }
    }
    if (nSlices == 0) { // 所有碱基都被过滤
      positions->returnAlignment(newAlignment);
      continue;
    }
    // 先登记切片数再入 LRU —— 否则 output 线程可能在计数前就消费掉某片。
    newAlignment->pendingPages.store(nSlices, std::memory_order_relaxed);
    cur = -1;
    int lo = 0;
    for (int i = 0; i < nb; ++i) {
      if (newAlignment->bases[i].remove)
        continue;
      int p =
          (int)((loc + newAlignment->bases[i].refPos - 1) / loadingBlockSize);
      if (p != cur) {
        if (cur != -1)
          positions->LRU.set(cur, PageSlice{newAlignment, lo, i});
        cur = p;
        lo = i;
      }
    }
    positions->LRU.set(cur, PageSlice{newAlignment, lo, nb});
  }
}
// 添加代码

// 一个渲染结果:整页文本 + 单调序号。done=true 表示某个 worker 已完成。
struct RenderItem {
  long long seq = 0;
  bool done = false;
  string text;
};

SafeQueue<RenderItem *> renderQueue;
const int renderQueueHighWater = 1024;

// 输出 worker:从 outputPool 取页节点,完成 A(页加载)/B(逐碱基累加)/C(整页渲染),
// 把 (seq, text) 交给 writer;取到哨兵则交付完成信号后退出。
void workerLoop(Positions *positions, string inputRefFileName) {
  PageWorker w;
  w.refFile.open(inputRefFileName, ios_base::in);

  DLinkedNode *node;
  while (true) {
    if (!positions->LRU.outputPool.popFront(node)) {
      this_thread::sleep_for(std::chrono::microseconds(1));
      continue;
    }
    if (node == nullptr) { // 哨兵:交付完成信号后退出
      RenderItem *done = new RenderItem();
      done->done = true;
      renderQueue.push(done);
      return;
    }
    const string &chr = node->vec.front().a->chromosome;
    positions->startload(w, chr, node->key);
    for (const PageSlice &slice : node->vec) {
      positions->appendPositions(w, slice, node->key);
      // 读段可能被多个页共享; 最后一片处理完才归还对象池。
      if (slice.a->pendingPages.fetch_sub(1, std::memory_order_acq_rel) == 1)
        positions->returnAlignment(slice.a);
    }
    positions->_renderPage(w);
    while (renderQueue.size() > renderQueueHighWater) { // 渲染侧背压
      this_thread::sleep_for(std::chrono::microseconds(10));
    }
    RenderItem *item = new RenderItem();
    item->seq = node->seq;
    item->text.swap(w.render); // 转移文本,避免拷贝
    renderQueue.push(item);
    positions->LRU.returnDLinkedNode(node);
  }
}

// writer:独占输出流,按 seq 升序落盘。乱序到达的结果先暂存,缺口补齐后按序写出。
void writerLoop(Positions *positions, string outputFileName, int nWorkers) {
  ostream *out_ = &cout;
  ofstream tableFile;
  if (!outputFileName.empty()) {
    tableFile.open(outputFileName, ios_base::out); // 覆盖模式,文件会被清空
    out_ = &tableFile;
    // 仅文件输出写表头,stdout 模式与旧版一致不写表头
    *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\t"
             "unconvertedBaseQualities\tunconvertedBaseCount\n";
  }

  long long nextToWrite = 0;
  map<long long, RenderItem *> pending;
  int doneCount = 0;
  while (doneCount < nWorkers) {
    RenderItem *item = nullptr;
    if (!renderQueue.popFront(item)) {
      this_thread::sleep_for(std::chrono::microseconds(1));
      continue;
    }
    if (item->done) {
      ++doneCount;
      delete item;
      continue;
    }
    if (item->seq == nextToWrite) {
      *out_ << item->text;
      delete item;
      ++nextToWrite;
      auto it = pending.find(nextToWrite);
      while (it != pending.end()) { // 冲刷已补齐的乱序暂存
        *out_ << it->second->text;
        delete it->second;
        pending.erase(it);
        ++nextToWrite;
        it = pending.find(nextToWrite);
      }
    } else {
      pending[item->seq] = item;
    }
  }
  tableFile.close();
}

int hisat_3n_table() {
  fprintf(stderr, "thread %s tid=%d\n", "main", gettid());
  positions = new Positions(nThreads, addedChrName,
                            removedChrName); // 给所有线程？分配一把锁？
  positions->LoadChromosomeNamesPos(refFileName);

  int k = nThreads < 1 ? 1 : nThreads; // worker 数 = -p 值;总线程 = k + main + writer
  vector<thread> workers;
  workers.reserve(k);
  for (int i = 0; i < k; i++)
    workers.emplace_back(workerLoop, positions, refFileName);
  thread writerThread(writerLoop, positions, outputFileName, k);

  ifstream inputFile;
  istream *alignmentFile = &cin;

  string *nextline = nullptr;
  // string *line;                // temporary string to get SAM line.

  long long int lastPos = 0; // the position on last SAM line. compare lastPos
                             // with samPos to make sure the SAM is sorted.

  if (!standardInMode) {
    inputFile.open(alignmentFileName, ios_base::in);
    alignmentFile = &inputFile;
  }
  forwardSam(alignmentFile, positions, nextline);
  while (readNewSamChromosome(alignmentFile, positions,
                              nextline)) { // read一整个染色体再统一处理

    processLinePool(positions); // 工作线程们   line=>LRU
    positions->LRU.popAllNodesFromTail();
    // 处理一群节点们
  }
  processLinePool(positions);
  positions->LRU.popAllNodesFromTail();

  if (!standardInMode) {
    inputFile.close();
  }
  // 等所有真实节点被 worker 取走,再每 worker 送一个哨兵
  while (!positions->LRU.outputPool.empty()) {
    this_thread::sleep_for(std::chrono::microseconds(100));
  }
  for (int i = 0; i < k; i++) {
    DLinkedNode *sentinel = nullptr;
    positions->LRU.outputPool.push(sentinel);
  }

  for (auto &t : workers)
    t.join();
  writerThread.join();
  delete positions;
  return 0;
}
int main(int argc, const char **argv) {
  int ret = 0;

  try {
    parseOptions(argc, argv);
    ret = hisat_3n_table();
  } catch (std::exception &e) {
    cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
    cerr << "Command: ";
    for (int i = 0; i < argc; i++)
      cerr << argv[i] << " ";
    cerr << endl;
    return 1;
  } catch (int e) {
    if (e != 0) {
      cerr << "Error: Encountered internal HISAT-3N exception (#" << e << ")"
           << endl;
      cerr << "Command: ";
      for (int i = 0; i < argc; i++)
        cerr << argv[i] << " ";
      cerr << endl;
    }
    return e;
  }

  return ret;
}
