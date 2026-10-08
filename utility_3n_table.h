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

#ifndef UTILITY_3N_TABLE_H
#define UTILITY_3N_TABLE_H
#include <algorithm>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <initializer_list>
#include <iosfwd>
#include <iostream>
#include <mutex>
#include <queue>
#include <unordered_map>
#include <vector>

using namespace std;

/**
 * return complement of input base.
 */
char asc2dnacomp[] = {
    /*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '-', 0, 0,
    /*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*  64 */ 0, 'T', 'V', 'G', 'H', 0, 0, 'C', 'D', 0, 0, 'M', 0, 'K', 'N', 0,
    /*    A   B   C   D           G   H           K       M   N */
    /*  80 */ 0, 0, 'Y', 'S', 'A', 0, 'B', 'W', 0, 'R', 0, 0, 0, 0, 0, 0,
    /*        R   S   T       V   W       Y */
    /*  96 */ 0, 'T', 'V', 'G', 'H', 0, 0, 'C', 'D', 0, 0, 'M', 0, 'K', 'N', 0,
    /*   a   b   c   d           g   h           k       m   n */
    /* 112 */ 0, 0, 'Y', 'S', 'A', 0, 'B', 'W', 0, 'R', 0, 0, 0, 0, 0, 0,
    /*        r   s   t       v   w       y */
    /* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/**
 * the simple data structure to bind quality score and position (on reference)
 * together.
 */
class PosQuality {
public:
  int readPos; // 0-based
  int refPos;  // 0-based
  char
      qual; // 碱基的测序质量值（Phred score），通常用 ASCII
            // 字符表示。它反映了测序仪对该碱基测序的可信程度：数字越大表示碱基越可靠
  bool converted;
  bool remove;

  PosQuality(int &inputPos) {
    readPos = inputPos;
    refPos = inputPos;
    remove = true;
  }

  void setQual(char &inputQual, bool inputConverted) {
    qual = inputQual;
    converted = inputConverted;
    remove = false;
  }
};

/**
 * the base class for string we need to search.
 */
class string_search {
public:
  int start;
  string s;
  int stringLen;

  void initialize() {
    start = 0;
    stringLen = 0;
    s.clear();
  }

  void loadString(string intputString) {
    s = intputString;
    stringLen = s.size();
    start = 0;
  }
};

/**
 * to store CIGAR string and search segments in it.
 */
class CIGAR : public string_search {
public:
  bool getNextSegment(int &len, char &symbol) {
    if (start == stringLen) {
      return false;
    }
    len = 0;
    int currentIndex = start;
    while (true) {
      if (isalpha(s[currentIndex])) {
        len = stoi(s.substr(start, currentIndex - start));
        symbol = s[currentIndex];
        start = currentIndex + 1;
        return true;
      }
      currentIndex++;
    }
  }
};

/**
 * to store MD tag and search segments in it.
 */
class MD_tag : public string_search {
public:
  bool getNextSegment(string &seg) {
    if (start >= stringLen) {
      return false;
    }
    seg.clear();
    int currentIndex = start;
    bool deletion = false;

    while (true) {
      if (currentIndex >= stringLen) {
        start = currentIndex + 1;
        return !seg.empty();
      }
      if (seg.empty() && s[currentIndex] == '0') {
        currentIndex++;
        continue;
      }
      if (isalpha(s[currentIndex])) {
        if (seg.empty()) {
          seg = s[currentIndex];
          start = currentIndex + 1;
          return true;
        } else {
          if (deletion) {
            seg += s[currentIndex];
            // currentIndex++;
          } else {
            start = currentIndex;
            return true;
          }
        }
      } else if (s[currentIndex] == '^') {
        if (seg.empty()) {
          seg = s[currentIndex];
          deletion = true;
        } else {
          start = currentIndex;
          return true;
        }
      } else { // number
        if (seg.empty()) {
          seg = s[currentIndex];
        } else {
          if (deletion || isalpha(seg.back())) {
            start = currentIndex;
            return true;
          } else {
            seg += s[currentIndex];
          }
        }
      }
      currentIndex++;
    }
  }
};

/**
 * simple safe queue
 */
template <typename T> class _SafeQueue {
private:
  mutex mutex_;
  queue<T> queue_;
  // condition_variable cv;
  bool closed = false;
  int queue_size = 0;

  string getReadName(string *line) {
    int startPosition = 0;
    int endPosition;

    endPosition = line->find("\t", startPosition);
    string readName = line->substr(startPosition, endPosition - startPosition);
    return readName;
  }

public:
  void pop() {
    std::unique_lock<std::mutex> lock(mutex_);
    queue_.pop();
  }

  T front() {
    std::unique_lock<std::mutex> lock(mutex_);
    T value = queue_.front();
    return value;
  }

  int size() {
    std::unique_lock<std::mutex> lock(mutex_);
    int s = queue_.size();
    return s;
  }

  /**
   * return true if the queue is not empty and pop front and get value.
   * return false if the queue is empty.
   */
  bool popFront(T &value) {
    std::unique_lock<std::mutex> lock(mutex_);
    bool isEmpty = queue_.empty();
    if (!isEmpty) {
      value = queue_.front();
      queue_.pop();
    }
    return !isEmpty;
  }

  void push(T value) {
    std::unique_lock<std::mutex> lock(mutex_);
    queue_.push(value);
  }

  bool empty() {
    std::unique_lock<std::mutex> lock(mutex_);
    bool check = queue_.empty();
    return check;
  }
};

/**
 * 线程安全队列。
 * 注:原实现为手写无锁队列(Michael-Scott),其 popFront 会 delete 被弹出的头节点,
 * 在多生产者/多消费者下存在 use-after-free;改为互斥锁实现以保证正确性。
 */
template <typename T> class SafeQueue {
  mutex mutex_;
  queue<T> queue_;

public:
  void push(T &value) {
    lock_guard<mutex> lock(mutex_);
    queue_.push(value);
  }

  bool popFront(T &result) {
    lock_guard<mutex> lock(mutex_);
    if (queue_.empty()) {
      return false;
    }
    result = queue_.front();
    queue_.pop();
    return true;
  }

  int size() {
    lock_guard<mutex> lock(mutex_);
    return (int)queue_.size();
  }

  bool empty() {
    lock_guard<mutex> lock(mutex_);
    return queue_.empty();
  }
};

/**
 * store one chromosome and it's stream position
 */
class ChromosomeFilePosition {
public:
  string chromosome;
  streampos linePos;
  int lineBases; // 每行碱基数(来自 .fai)
  int lineWidth; // 每行字节数,含换行(来自 .fai)
  // streampos endlinePos;
  ChromosomeFilePosition(string inputChromosome, streampos inputstartPos,
                         int inputLineBases, int inputLineWidth) {
    chromosome = inputChromosome;
    linePos = inputstartPos;
    lineBases = inputLineBases;
    lineWidth = inputLineWidth;
    // endlinePos = inputendPos;
  }

  bool operator<(const ChromosomeFilePosition &in) const {
    return chromosome < in.chromosome;
  }
};

/**
 * store all chromosome and it's stream position
 */
class ChromosomeFilePositions {
public:
  vector<ChromosomeFilePosition> pos;

  /**
   * input the chromosome name and it's streamPos, if it is not in pos, add it.
   */
  void append(string &chromosome, streampos &linePos, int lineBases,
              int lineWidth) {
    pos.push_back(
        ChromosomeFilePosition(chromosome, linePos, lineBases, lineWidth));
  }

  /**
   * make binary search on pos for target chromosome name
   */
  int findChromosome(const string &targetChromosome, int start, int end) const {
    if (start <= end) {
      int middle = (start + end) / 2;
      if (pos[middle].chromosome == targetChromosome) {
        return middle;
      }
      if (pos[middle].chromosome > targetChromosome) {
        return findChromosome(targetChromosome, start, middle - 1);
      }
      return findChromosome(targetChromosome, middle + 1, end);
    } else {
      // cannot find the chromosome! throw!
      cerr << "Cannot find the chromosome: " << targetChromosome
           << " in reference file." << endl;
      throw 1;
    }
  }

  /**
   * 返回匹配到的条目指针(供 per-worker 缓存复用),未命中则抛异常。
   */
  const ChromosomeFilePosition *findEntry(const string &targetChromosome) const {
    int index = findChromosome(targetChromosome, 0, (int)pos.size() - 1);
    return &pos[index];
  }

  /**
   * given targetChromosome name, return its streampos
   */
  streampos getChromosomePosInRefFile(const string &targetChromosome) {
    int index = findChromosome(targetChromosome, 0, pos.size() - 1);
    assert(pos[index].chromosome == targetChromosome);
    return pos[index].linePos;
  }

  /**
   * byte = linePos + (basePos / lineBases) * lineWidth + (basePos % lineBases)
   */
  streampos byteOffsetOf(const string &targetChromosome, long long int basePos) {
    int index = findChromosome(targetChromosome, 0, pos.size() - 1);
    const ChromosomeFilePosition &c = pos[index];
    streamoff within = (streamoff)((basePos / c.lineBases) * c.lineWidth +
                                   (basePos % c.lineBases));
    return c.linePos + within;
  }

  int lineBasesOf(const string &targetChromosome) {
    int index = findChromosome(targetChromosome, 0, pos.size() - 1);
    return pos[index].lineBases;
  }

  /**
   * sort the pos by chromosome name
   */
  void sort() { std::sort(pos.begin(), pos.end()); }
};
#pragma one
class Alignment;
/**
 * 一个 (页, 读段) 切片: 该读段落在本页的碱基下标区间 [lo, hi)。
 * 读段对象由它跨越的所有页共享, 切片避免了为每页复制一份读段。
 */
struct PageSlice {
  Alignment *a;
  int lo;
  int hi;
};
struct DLinkedNode {
  int key;
  long long seq = 0;     // 输出序号:入 outputPool 前由 main 单生产者单调分配
  vector<PageSlice> vec; // 本页上的读段切片
  DLinkedNode *prev;
  DLinkedNode *next;
  DLinkedNode(int k = 0)
      : key(k), prev(nullptr), next(nullptr) {
  } // std::vector 自带默认构造函数，不用手动初始化
  void initialize() {
    vector<PageSlice>().swap(vec);
    key = 0;
    seq = 0;
    prev = nullptr;
    next = nullptr;
  }
};

class LRUCache {
private:
  unordered_map<int, DLinkedNode *> cache;
  int capacity;
  int count;
  DLinkedNode *head;
  DLinkedNode *tail;
  mutable mutex mtx; // 🔒全局互斥锁
  long long nextSeq = 0; // 输出序号计数器(仅在 main 持 mtx 时递增)

public:
  SafeQueue<DLinkedNode *> outputPool;
  SafeQueue<DLinkedNode *> freeDLinkedNodePool;

  LRUCache(int capacity) : capacity(capacity), count(0) {
    head = new DLinkedNode();
    tail = new DLinkedNode();
    head->next = tail;
    tail->prev = head;
  }
  ~LRUCache() {
    DLinkedNode *node;
    while (freeDLinkedNodePool.popFront(node)) {
      delete node;
    }

    while (outputPool.popFront(node)) {
      delete node;
    }
    delete head;
    delete tail;
  }
  void set(int key, const PageSlice &slice) {
    lock_guard<mutex> lock(mtx);
    if (cache.find(key) == cache.end()) {
      DLinkedNode *newNode;
      getFreeDLinkedNode(newNode, key); // safequeue
      cache[key] = newNode;
      newNode->vec.emplace_back(slice);
      addNode(newNode);
      ++count;

      if (count > capacity) {
        DLinkedNode *tailNode = popTail();
        cache.erase(tailNode->key);
        // delete tailNode;
        --count;
        tailNode->seq = nextSeq++; // 入队前编号,复现单线程消费顺序
        outputPool.push(tailNode);
        return; // TODO:有弹出，需要输出,接着回收
      }
    } else {
      DLinkedNode *node = cache[key];
      node->vec.emplace_back(slice);
      moveToHead(node);
    }
    return;
  }
  void popAllNodesFromTail() {
    lock_guard<mutex> lock(mtx); // 加锁，防止并发访问

    DLinkedNode *node = tail->prev;
    while (node != head) { // 从尾部往前遍历直到head
      DLinkedNode *prev = node->prev;
      removeNode(node);
      cache.erase(node->key);
      --count;
      node->seq = nextSeq++; // 入队前编号,复现单线程消费顺序
      outputPool.push(node);
      node = prev;
    }
  }

  void getFreeDLinkedNode(DLinkedNode *&newDLinkedNode,
                          int key) { // 对一个指向 Alignment
                                     // 的指针的引用,只用*无法改变外部
    if (freeDLinkedNodePool.popFront(newDLinkedNode)) {
      newDLinkedNode->key = key;
      return;
    } else {
      newDLinkedNode = new DLinkedNode(key);
    }
  }
  void returnDLinkedNode(DLinkedNode *newDLinkedNode) {
    newDLinkedNode->initialize();
    //: 写alignment的clear函数,已经写好了
    freeDLinkedNodePool.push(newDLinkedNode);
  }

private:
  void addNode(DLinkedNode *node) {
    node->prev = head;
    node->next = head->next;
    head->next->prev = node;
    head->next = node;
  }

  void removeNode(DLinkedNode *node) {
    DLinkedNode *prev = node->prev;
    DLinkedNode *next = node->next;
    prev->next = next;
    next->prev = prev;
  }

  void moveToHead(DLinkedNode *node) {
    removeNode(node);
    addNode(node);
  }

  DLinkedNode *popTail() {
    DLinkedNode *res = tail->prev;
    removeNode(res);
    return res;
  }
};

#endif // UTILITY_3N_TABLE_H
