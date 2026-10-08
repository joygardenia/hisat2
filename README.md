HISAT-3N
============

> **本仓库是 HISAT-3N 的性能优化分支**
> 将 `hisat-3n-table` 从「全基因组逐碱基扫描」重构为「按页按需加载 + 多 worker Pipeline + 保序输出」
> 冷缓存口径加速 **132×**(首次运行,含构建索引)/ **766×**(索引已就绪);
> 高密度数据(10% 覆盖)下原始 vs 优化加速 **29.8× / 28.4×**(4 / 8 线程)
> 上游项目:https://github.com/DaehwanKimLab/hisat2 (GPLv3)

## 本分支做了什么

`hisat-3n-table` 从排序后的 SAM 与参考基因组,计算每个位点的碱基转化统计。

**原始实现的瓶颈不在 IO,在 CPU。** 主循环随着 SAM 位置推进,顺序扫完整份参考基因组:

```cpp
while (samPos > reloadPos) {   // SAM 跨越全部染色体,reloadPos 被推着走完整个基因组
    positions->loadMore();
    reloadPos += loadingBlockSize;
}
```

`loadMore()` 对扫过的**每一个碱基**执行:

- `toupper()` 逐字符大写
- `getFreePosition()` 取对象 + `refPositions.push_back()` 容器增长

主线程还要用 `sleep_for(1μs)` spin-wait 等 worker。因此原版成本是 **O(参考基因组大小)**,与 SAM 比对条数无关 —— 本例中仅 663 条比对,却要跑 146 秒。

优化后的做法:

| 优化 | 做法 | 效果 |
|---|---|---|
| 按页按需加载 | 按 SAM 命中的 ref 页号直接 `seekg`,取代顺序推进 | 复杂度 O(基因组) → O(命中范围) |
| 预分配数组 | `ref + (index_group++)` 取代每碱基 `new` + `push_back` | 热路径零堆分配 |
| 多 worker Pipeline | main 解析/分页,k 个 worker 并行 Ref-Load/Compute/Output(页独占) | 输出阶段 4~8 路并行 |
| 保序 writer | 页节点带单调序号,单 writer 按序落盘(乱序暂存) | 多线程输出与原单线程逐字节一致 |
| LRU 页缓存 | 以 ref 页号为 key 聚合 read,按染色体边界批量淘汰 | 同页 read 合并 |
| 三级对象池 | `string` / `Alignment` / `DLinkedNode` 复用 | 消除跨线程分配开销 |
| 线程安全队列 | 阶段之间以互斥锁队列传递(修正原无锁队列的 use-after-free) | 取代 spin-wait、保证正确性 |
| 去掉位点锁 | 页独占后删除每碱基 `lock/unlock` | 加速逐碱基累加 |
| 染色体偏移索引 | 扫描结果落盘为 samtools 兼容的 `<ref>.fai`,取代每次运行的全量扫描 | 稳态下 ref 读取 3 GB → 6.4 KB |
| 有界内存扫描 | 索引构建改为 8 MB 分窗 + `MADV_DONTNEED` | 构建期 maxRSS 3083 MB → 14 MB |

### 性能

数据集:`SRR26378479_filter.dump`(663 条比对)+ GRCh38(3.0 GB)。
每轮运行前 `sync` + `drop_caches` 冷启动,取 3 次运行的中位数。

| 口径 | 原始 | 优化后 | 加速比 | maxRSS |
|---|---|---|---|---|
| 首次运行(含构建 `<ref>.fai` 索引) | 145.62 s | **1.10 s** | **132×** | 401 → 14 MB |
| 索引已就绪 | 145.62 s | **0.19 s** | **766×** | 401 → 6 MB |

**索引成本的计算口径**:原始实现每次运行都必须顺序扫完整份参考基因组,才能确定各染色体的字节偏移;
优化后把这一步的产物落盘(`.fai`),首次构建之后即可复用。同一份参考基因组跑 N 次的总成本:

    原始:   N × 145.62 s
    优化后: 1.10 s + (N − 1) × 0.19 s

N=1 → 132×,N=3 → 295×,N 增大时渐近 **766×**。

**为什么原版是 CPU 密集型而非 IO 密集型**:原版实测 CPU 占用 108%,全程跑满;
瓶颈是 31 亿次 `toupper` + 每碱基构造 `Position`,外加一千万次 spin-wait 上下文切换。
相比之下,读 3 GB 参考基因组按 NVMe 带宽只需约 1 秒 —— 磁盘从来不是瓶颈。

**正确性**:两版本输出均为 3122 行,排序后逐字节一致。

### 多线程性能(高密度)

数据集:GRCh38 上合成的 **10% 覆盖**数据(stride=1000,**2,948,601 条比对**,输出 6,011 万行)。
口径:governor=performance,每轮 `sync` + `drop_caches` 冷启动,取 2 次中位数。

| 线程 | 原始 | 优化后 | 加速比 | maxRSS |
|---|---|---|---|---|
| `-p 4` | 370.58 s | **12.43 s** | **29.8×** | 404 → 306 MB |
| `-p 8` | 388.71 s | **13.70 s** | **28.4×** | 403 → 296 MB |

- 输出经 `cmp` 校验,原始版与优化版**逐字节一致**。
- 原版加线程**更慢**(4→8:370.6→388.7 s);优化版 4→8 基本持平(12.4→13.7 s),
  已达 main 线程(~13.6 s)天花板 —— **4~8 worker 为甜点**。
- 优化版为**稳态**口径(保留已有 `.fai`,不重建索引)。

### 测试环境

| 项 | 值 |
|---|---|
| CPU | Intel Core i5-13500H (16 threads),`governor=performance` |
| 内存 | 15 GiB |
| 磁盘 | NVMe SSD |
| OS / 编译器 | Manjaro Linux 6.12.103 / GCC 16.1.1 |
| 测量方法 | 冷口径每轮运行前 `sync` + `echo 3 > drop_caches`;取中位数(N=3) |

### 已知代价与局限

- **首次运行需要构建索引**:约 1.10 s(构建 + 运行),之后每次 0.19 s。
  索引为 samtools 兼容格式,可直接复用 `samtools faidx <ref.fa>` 的产物,无需本程序生成。
- **`-p/--threads` 生效**:映射为 output worker 数(缺省 1;总线程 = worker + main + writer)。
  4~8 worker 为甜点,超过后受 main 线程(~13.6 s)与内存带宽限制,收益递减。
  原版 `-p` 相反:其瓶颈在主线程,增加线程只会加剧锁竞争(实测原版 4→8 线程更慢)。

### 构建与运行

    git clone https://github.com/joygardenia/hisat2.git
    cd hisat2
    git checkout release
    make hisat-3n-table

    ./hisat-3n-table -m -p 8 \
      --alignments <sorted.sam> \
      --ref <reference.fa> \
      --output-name /dev/stdout \
      --base-change C,T

首次运行会在参考基因组旁自动生成 `<reference.fa>.fai`(6.4 KB),之后运行直接复用。
如需把索引构建移出计时窗口,可先用 `samtools faidx <reference.fa>` 预建(格式兼容)。

### 改动文件

| 文件 | 职责 |
|---|---|
| `hisat_3n_table.cpp` | `hisat-3n-table` 入口:N worker + 1 保序 writer 调度、`-p` 映射、`.fai` 初始化 |
| `position_3n_table.h` | `Positions` / `PageWorker`:染色体偏移索引(`.fai` 读写)、per-worker 按页加载 ref、逐位点统计 |
| `utility_3n_table.h` | `SafeQueue` / `LRUCache`(页节点序号) / `ChromosomeFilePositions` / 三级对象池 |
| `alignment_3n_table.h` | SAM 记录解析、`front_page` 页号计算 |

---

<!-- 以下为上游 HISAT-3N 原始文档,未作修改 -->

Overview
-----------------
HISAT-3N (hierarchical indexing for spliced alignment of transcripts - 3 nucleotides)
is an ultrafast and memory-efficient sequence aligner designed for nucleotide conversion
sequencing technologies. HISAT-3N index contains two HISAT2 indexes which require memory small:
for the human genome, it requires 9 GB for standard 3N-index and 10.5 GB for repeat 3N-index.
The repeat 3N-index could be used to align one read to thousands position 3 times faster standard 3N-index.
HISAT-3N is developed based on [HISAT2],
which is particularly optimized for RNA sequencing technology. HISAT-3N support both strand-specific and non-strand reads.
HISAT-3N can be used for any base-converted sequencing reads include [BS-seq], [SLAM-seq], [scBS-seq], [scSLAM-seq], and [TAPS].
See the [HISAT-3N] website for more information.

[HISAT2]:https://github.com/DaehwanKimLab/hisat2
[BS-seq]: https://en.wikipedia.org/wiki/Bisulfite_sequencing
[SLAM-seq]: https://www.nature.com/articles/nmeth.4435
[scBS-seq]: https://www.nature.com/articles/nmeth.3035
[scSLAM-seq]: https://www.nature.com/articles/s41586-019-1369-y
[TAPS]: https://www.nature.com/articles/s41587-019-0041-2
[HISAT-3N]:https://daehwankimlab.github.io/hisat2/hisat-3n


Getting started
============
HISAT-3N requires a 64-bit computer running either Linux or Mac OS X and at least 16 GB of RAM.

A few notes:

1. Building the standard 3N index requires 16GB of RAM or less.
2. Building the repeat 3N index requires 256GB of RAM.
3. The alignment process using either the standard or repeat index requires less than 16GB of RAM.
4. [SAMtools] is required to sort SAM files in order to generate a HISAT-3N table.

Install
------------

    git clone https://github.com/DaehwanKimLab/hisat2.git hisat-3n
    cd hisat-3n
    git checkout -b hisat-3n origin/hisat-3n
    make

Build a HISAT-3N index with `hisat-3n-build`
-----------
`hisat-3n-build` builds a 3N-index, which contains two hisat2 indexes, from a set of DNA sequences. For standard 3N-index,
each index contains 16 files with suffix `.3n.*.*.ht2`.
For repeat 3N-index, there are 16 more files in addition to the standard 3N-index, and they have the suffix
`.3n.*.rep.*.ht2`.
These files constitute the hisat-3n index and no other file is needed to alignment reads to the reference.

* `--base-change <chr1,chr2>` argument is required for `hisat-3n-build` and `hisat-3n`.   
  Provide which base is converted in the sequencing process to another base. Please enter
  2 letters separated by ',' for this argument. The first letter(chr1) should be the converted base, the second letter(chr2) should be
  the converted to base. For example, during slam-seq, some 'T' is converted to 'C',
  please enter `--base-change T,C`. During bisulfite-seq, some 'C' is converted to 'T', please enter `--base-change C,T`.
* Different conversion types may build the same hisat-3n index. Please check the table below for more detail.
  Once you build the hisat-3n index with C to T conversion (for example BS-seq).
  You can align the T to C conversion reads (for example SLAM-seq reads) with the same index.


| Conversion Types                   | HISAT-3N index suffix         |
  |:----------------------------------:|:-----------------------------:|
|C -> T<br>T -> C<br>A -> G<br>G -> A|.3n.CT.\*.ht2 <br>.3n.GA.\*.ht2|
|A -> C<br>C -> A<br>G -> T<br>T -> G|.3n.AC.\*.ht2 <br>.3n.TG.\*.ht2|
|A -> T<br>T -> A                    |.3n.AT.\*.ht2 <br>.3n.TA.\*.ht2|
|C -> G<br>G -> C                    |.3n.CG.\*.ht2 <br>.3n.GC.\*.ht2|

#### Examples:
    # Build the standard HISAT-3N index (with C to T conversion):  
    hisat-3n-build --base-change C,T genome.fa genome
  
    # Build the repeat HISAT-3N index (with T to C conversion, require 256 GB memory for human genome index):  
    hisat-3n-build --base-change T,C --repeat-index genome.fa genome

It is optional to make the graph index and add SNP or spice site information to the index, to increase the alignment accuracy.
The graph index building may require more memory than the linear index building.
For more detail, please check the [HISAT2 manual].

[HISAT2 manual]:https://daehwankimlab.github.io/hisat2/manual/

#### Examples:
    # Build the standard HISAT-3N index integrated index with SNP information
    hisat-3n-build --base-change C,T --snp genome.snp genome.fa genome 
    
    # Build the standard HISAT-3N integrated index with splice site information
    hisat-3n-build --base-change C,T --ss genome.ss --exon genome.exon genome.fa genome 
    
    # Build the repeat HISAT-3N index integrated index with SNP information
    hisat-3n-build --base-change C,T --repeat-index --snp genome.snp genome.fa genome 
    
    # Build the repeat HISAT-3N integrated index with splice site information
    hisat-3n-build --base-change C,T --repeat-index --ss genome.ss --exon genome.exon genome.fa genome 


Alignment with `hisat-3n`
------------
After building the HISAT-3N index, you are ready to use `hisat-3n` for alignment.
HISAT-3N has the same set of parameters as in HISAT2 with some additional arguments. Please refer to the [HISAT2 manual] for more details.

For the human reference genome, HISAT-3N requires about 9GB for alignment with the standard 3N-index and 10.5GB for the repeat 3N-index.

* `--base-change <nt1,nt2>`  
  Specify the nucleotide conversion type (e.g., C to T in bisulfite-sequencing reads). The parameter option is two characters separated by ','.  Type the original nucleotide for the first character (nt1) and type the converted nucleotide as the second character (nt2). For example, if performing [SLAM-seq] where some 'T's are converted to 'C's, input `--base-change T,C`.
  As another example, if performing bisulfite-seq, where some 'C's are converted to 'T's, please input `--base-change C,T`.
  If you want to align non-converted reads to the regular HISAT2 index, then omit this command.

* `--index/-x <hisat-3n-idx>`  
  Specify the index file basename for HISAT-3N.  The basename is the name of the index files up to but not including the suffix `.3n.*.*.ht2` / etc.
  For example, if you build your index with basename 'genome' using a HISAT-3N-build, please input `--index genome`.

* `--directional-mapping`  
  Make directional mapping. Please use this option only if your sequencing reads are generated from a strand-specific library. 
  The directional mapping mode is about 2x faster than the standard (non-directional) mapping mode.

* `--repeat-limit <int>`  
  You can set up the number of alignments to be checked for each repeat alignment. You may increase the number to direct hisat-3n
  to output more, if a read has multiple mapping locations. We suggest that you limit the repeat number for paired-end read alignment to no more
  than 1,000,000. default: 1000.

* `--unique-only`  
  Only output uniquely aligned reads.


#### Examples:
* Single-end [SLAM-seq] read (T to C conversion) alignment with standard 3N-index:  
  `hisat-3n --index genome -f -U read.fa -S output.sam --base-change T,C`

* Paired-end strand-specific bisulfite-seq read (C to T conversion) alignment with repeat 3N-index:   
  `hisat-3n --index genome -f -1 read_1.fa -2 read_2.fa -S output.sam --base-change C,T --directional-mapping`

* Single-end TAPS reads (C to T conversion) alignment with repeat 3N-index and only output unique aligned results:   
  `hisat-3n --index genome -q -U read.fq -S output.sam --base-change C,T --unique`



#### Extra SAM tags generated by HISAT-3N:

* `Yf:i:<N>`: Number of conversions detected in the read.
* `Zf:i:<N>`: Number of un-converted bases are detected in the read. Yf + Zf = total number of bases which can be converted in the read sequence.
* `YZ:A:<A>`: The value `+` or `–` indicates the read is mapped to REF-3N (`+`) or REF-RC-3N (`-`), respectively.

Generate a 3N-conversion-table with `hisat-3n-table`
------------
### Preparation

To generate a 3N-conversion-table, users need to sort the `hisat-3n` generated SAM alignment file.

[SAMtools] is required for this sorting process.

Use `samtools sort` to convert the SAM file into a sorted SAM file.

    samtools sort output.sam -o output_sorted.sam -O sam

Generate 3N-conversion-table with `hisat-3n-table`:

### Usage
    hisat-3n-table [options]* --alignments <alignmentFile> --ref <refFile> --base-change <char1,char2>

#### Main arguments
* `--alignments <alignmentFile>`   
  SORTED SAM file. Please enter `-` for standard input.

* `--ref <refFile>`  
  The reference genome file (FASTA format) for generating HISAT-3N index.

* `--output-name <outputFile>`  
  Filename to write 3N-conversion-table (tsv format) to.  By default, table is written to the “standard out” or “stdout” filehandle (i.e. the console).

* `--base-change <char1,char2>`  
  The base-change rule. User should enter the exact same `--base-change` arguments in hisat-3n.
  For example, please enter `--base-change C,T` for bisulfite sequencing reads.

#### Input options
* `-u/--unique-only`  
  Only count the unique aligned reads into 3N-conversion-table.

* `-m/--multiple-only`  
  Only count the multiple aligned reads into 3N-conversion-table.

* `-c/--CG-only`  
  Only count the CpG sites in reference genome. This option is designed for bisulfite sequencing reads.

* `--added-chrname`  
  Please add this option if you use `--add-chrname` during `hisat-3n` alignment.
  During `hisat-3n` alignment, the prefix "chr" is added in front of chromosome name and shows on SAM output, when user choose `--add-chrname`.
  `hisat-3n-table` cannot find the chromosome name on reference because it has an additional "chr" prefix. This option is to help `hisat-3n-table`
  find the matching chromosome name on reference file. The 3n-table provides the same chromosome name as SAM file.

* `--removed-chrname`  
  Please add this option if you use `--remove-chrname` during `hisat-3n` alignment.
  During `hisat-3n` alignment, the prefix "chr" is removed in front of chromosome name and shows on SAM output, when user choose `--remove-chrname`.
  `hisat-3n-table` cannot find the chromosome name on reference because it has no "chr" prefix. This option is to help `hisat-3n-table`
  find the matching chromosome name on reference file. The 3n-table provides the same chromosome name as SAM file.

#### Other options:
* `-p/--threads <int>`  
  Launch `int` parallel threads (default: 1) for table building.

* `-h/--help`  
  Print usage information and quit.

#### Examples:
    # Generate the 3N-conversion-table for bisulfite sequencing data:  
      hisat-3n-table -p 16 --alignments sorted_alignment_result.sam --ref genome.fa --output-name output.tsv --base-change C,T
    
    # Generate the 3N-conversion-table for TAPS data and only count base in CpG site and uniquely aligned:  
      hisat-3n-table -p 16 --alignments sorted_alignment_result.sam --ref genome.fa --output-name output.tsv --base-change C,T --CG-only --unique-only
    
    # Generate the 3N-conversion-table for bisulfite sequencing data from sorted BAM file:  
      samtools view -h sorted_alignment_result.bam | hisat-3n-table --ref genome.fa --alignments - --output-name output.tsv --base-change C,T
    
    # Generate the 3N-conversion-table for bisulfite sequencing data from unsorted BAM file:  
      samtools sort alignment_result.bam -O sam | hisat-3n-table --ref genome.fa --alignments - --output-name output.tsv --base-change C,T


#### Note:
There are 7 columns in the 3N-conversion-table:

1. `ref`: the chromosome name.
2. `pos`: 1-based position in `ref`.
3. `strand`: '+' for forward strand. '-' for reverse strand.
4. `convertedBaseQualities`: the qualities of the converted bases in read-level measurement. The length of this string is equal to the number of converted bases.
5. `convertedBaseCount`: the number of distinct read positions where converted bases in read-level measurements were found.
   this number is equal to the length of convertedBaseQualities.
6. `unconvertedBaseQualities`: the qualities of the unconverted bases in read-level measurement. The length of this string is equal to the number of unconverted bases in read-level measurement.
7. `unconvertedBaseCount`: the number of distinct read positions where unconverted bases in read-level measurements were found.
   this number is equal to the length of unconvertedBaseQualities.

##### Sample 3N-conversion-table:
    ref    pos    strand    convertedBaseQualities    convertedBaseCount    unconvertedBaseQualities    unconvertedBaseCount
    1      11874  +         FFFFFB<BF<F               11                                                0
    1      11877  -         FFFFFF<                   7                                                 0
    1      11878  +         FFFBB//F/BB               11                                                0
    1      11879  +                                   0                     FFFBB//FB/                  10
    1      11880  -         F                         1                     FFFF/                       5
[SAMtools]:        http://samtools.sourceforge.net

Publication
============

* HISAT-3N   
  Zhang, Y., Park, C., Bennett, C., Thornton, M. and Kim, D. [Rapid and accurate alignment of nucleotide conversion sequencing reads with HISAT-3N](https://doi.org/10.1101/gr.275193.120). _Genome Research_ **31(7)**: 1290-1295 (2021)


* HIAST2   
  Kim, D., Paggi, J.M., Park, C. _et al._ [Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype](https://doi.org/10.1038/s41587-019-0201-4). _Nat Biotechnol_ **37**, 907–915 (2019)  
