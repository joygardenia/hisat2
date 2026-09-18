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

template <typename T> class SafeQueue {
  struct Node {
    T data;
    std::atomic<Node *> next;
    Node(const T &val, Node *n) : data(val) {
      next.store(n, std::memory_order_relaxed);
    }
    Node() : next(nullptr) {} // 用于 dummy 节点
  };
  std::atomic<Node *> head, tail;
  std::atomic<int> queue_size{0};

public:
  SafeQueue() {
    Node *dummy = new Node;
    head = tail = dummy;
  }

  void push(T &value) {
    Node *node = new Node(value, nullptr);
    Node *old_tail;
    while (true) {
      old_tail = tail.load(std::memory_order_acquire);
      Node *next = old_tail->next.load(std::memory_order_acquire);
      if (next == nullptr) {
        if (old_tail->next.compare_exchange_weak(next, node))
          break;
      } else {
        tail.compare_exchange_weak(old_tail, next);
      }
    }
    tail.compare_exchange_weak(old_tail, node);
    queue_size.fetch_add(1, std::memory_order_relaxed);
  }

  bool popFront(T &result) {
    Node *old_head;
    while (true) {
      old_head = head.load(std::memory_order_acquire);
      Node *next = old_head->next.load(std::memory_order_acquire);
      if (next == nullptr) {
        // queue_size.fetch_sub(1);   队列为空时错误地减少
        // size,没有弹出元素，想当然了
        return false;
      }
      if (head.compare_exchange_weak(old_head, next)) {
        result = next->data;
        delete old_head;
        queue_size.fetch_sub(1, std::memory_order_relaxed);
        return true;
      }
    }
  }

  int size() { return queue_size.load(std::memory_order_acquire); }

  bool
  empty() { // return head.load() == tail.load();逻辑对，但在无锁情况下，这两个
            // load 没同步顺序，可能被 CPU 乱序
    return head.load(std::memory_order_acquire) ==
           tail.load(std::memory_order_acquire);
  }
};

/**
 * store one chromosome and it's stream position
 */
class ChromosomeFilePosition {
public:
  string chromosome;
  streampos linePos;
  // streampos endlinePos;
  ChromosomeFilePosition(string inputChromosome, streampos inputstartPos) {
    chromosome = inputChromosome;
    linePos = inputstartPos;
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
  void append(string &chromosome, streampos &linePos) {
    pos.push_back(ChromosomeFilePosition(chromosome, linePos));
  }

  /**
   * make binary search on pos for target chromosome name
   */
  int findChromosome(string &targetChromosome, int start, int end) {
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
   * given targetChromosome name, return its streampos
   */
  streampos getChromosomePosInRefFile(string &targetChromosome) {
    int index = findChromosome(targetChromosome, 0, pos.size() - 1);
    assert(pos[index].chromosome == targetChromosome);
    return pos[index].linePos;
  }

  /**
   * sort the pos by chromosome name
   */
  void sort() { std::sort(pos.begin(), pos.end()); }
};
#pragma one
class Alignment;
struct DLinkedNode {
  int key;
  // int value; // 此处用refpositions的类
  vector<Alignment *> vec; // no 这是一个vector,存放alingment
  DLinkedNode *prev;
  DLinkedNode *next;
  DLinkedNode(int k = 0)
      : key(k), prev(nullptr), next(nullptr) {
  } // std::vector 自带默认构造函数，不用手动初始化
  void initialize() {
    vector<Alignment *>().swap(vec);
    key = 0;
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
  void set(int key, Alignment *alignment) {
    lock_guard<mutex> lock(mtx);
    if (cache.find(key) == cache.end()) {
      DLinkedNode *newNode;
      getFreeDLinkedNode(newNode, key); // safequeue
      cache[key] = newNode;
      newNode->vec.emplace_back(alignment);
      addNode(newNode);
      ++count;

      if (count > capacity) {
        DLinkedNode *tailNode = popTail();
        cache.erase(tailNode->key);
        // delete tailNode;
        --count;
        outputPool.push(tailNode);
        return; // TODO:有弹出，需要输出,接着回收
      }
    } else {
      DLinkedNode *node = cache[key];
      node->vec.emplace_back(alignment);
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
      outputPool.push(node);
      node = prev;
    }
  }

  void getFreeDLinkedNode(DLinkedNode *&newDLinkedNode,
                          int key) { // 对一个指向 Alignment
                                     // 的指针的引用,只用*无法改变外部
    if (freeDLinkedNodePool.popFront(newDLinkedNode)) {
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
