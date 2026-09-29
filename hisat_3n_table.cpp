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
int skipSize = 61;

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
void processLinePool(Positions *positions) {
  Alignment *newAlignment;
  while (positions->linePool.popFront(line)) {
    positions->getFreeAlignment(newAlignment);
    bool parsed = newAlignment->parse(line); // 解析
    positions->returnLine(line);
    if (!parsed) {
      positions->returnAlignment(newAlignment);
      continue;
    }
    positions->LRU.set(newAlignment->front_page, newAlignment);
    if (newAlignment->middleRefPos != -1) {
      Alignment *newAlignment_2;
      positions->getFreeAlignment(newAlignment_2);
      *newAlignment_2 = *newAlignment;
      newAlignment_2->isFirstPage = false;
      newAlignment_2->front_page += 1;
      positions->LRU.set(newAlignment_2->front_page, newAlignment_2);
    }
  }
}
// 添加代码
bool working = true;
void processSingleNode(Positions *positions,
                       string inputRefFileName) { // 输出线程
  fprintf(stderr, "thread %s tid=%d\n", "output", gettid());
  positions->refFile.open(inputRefFileName, ios_base::in);
  positions->LoadChromosomeNamesPos(inputRefFileName);

  DLinkedNode *node;
  while (working) { // 处理多个节点
    if (positions->LRU.outputPool.popFront(node)) {
      positions->loadTestChromosome(node->vec.front()->chromosome);
      positions->startload(
          node->key); // FIX 加载一页碱基 ,单条染色体全部加载进去了，此处会oom
      for (Alignment *newAlignment : node->vec) {
        positions->appendPositions(*newAlignment);
        positions->returnAlignment(newAlignment);
      }
      positions->_moveAllToprint(outputFileName);
    }
  }
}

int hisat_3n_table() {
  fprintf(stderr, "thread %s tid=%d\n", "main", gettid());
  positions = new Positions(nThreads, addedChrName,
                            removedChrName); // 给所有线程？分配一把锁？
  // main function, initially 2 load loadingBlockSize (2,000,000) bp of
  // reference, set reloadPos to 1 loadingBlockSize, then load SAM data. when
  // the samPos larger than the reloadPos load 1 loadingBlockSize bp of
  // reference. when the samChromosome is different to current chromosome,
  // finish all sam position and output all.
  thread outputThread;
  outputThread = thread(processSingleNode, positions, refFileName);
  ifstream inputFile;
  istream *alignmentFile = &cin;

  string *nextline = nullptr;
  // string *line;                // temporary string to get SAM line.

  long long int reloadPos; // the position in reference that we need to
                           // reload.反映参考基因组中的物理位置(基因组坐标级别）
  long long int lastPos = 0; // the position on last SAM line. compare lastPos
                             // with samPos to make sure the SAM is sorted.

  if (!standardInMode) {
    inputFile.open(alignmentFileName, ios_base::in);
    alignmentFile = &inputFile;
  }
  // 添加代码
  positions->writeHeader(outputFileName);
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
  // 收尾
  while (!positions->LRU.outputPool.empty()) {
    this_thread::sleep_for(std::chrono::microseconds(100));
  }
  working = false;
  outputThread.join();
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
