#!/usr/bin/env python3
"""跨多页读段(RNA-seq 大内含子)越界崩溃的最小复现。

页大小 = loadingBlockSize = 6000 bp。
read1 的 20000 bp 内含子使其映射横跨第 0 页与第 3 页 —— 旧实现只把读段切成
2 份(假设最多跨 2 页),第二份被错误地摆到第 1 页,于是
`ref + (index + refPos)` 越过 6000 槽缓冲 -> 野指针 -> mutex::lock 崩溃。

生成到 <outdir>:
  ref.fa            chr1, 30000 bp, "ACGT" 重复, 行宽 60
  case.sam          3 条读段(见下)
  expected_pos.list 独立 oracle: 参考中 C 碱基在映射区间内的 1-based 位置

ref.fa.fai 由 samtools faidx 生成(见 run.sh)。
"""
import os
import sys

CHR = "chr1"
LEN = 30000
LINE = 60


def ref_seq():
    pat = "ACGT"
    return (pat * (LEN // len(pat) + 1))[:LEN]


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(
        os.path.abspath(__file__))
    seq = ref_seq()

    # (name, 0-based mapped windows, cigar)
    reads = [
        ("read1", [(0, 50), (20050, 20100)], "50M20000N50M"),     # 跨页 0 -> 3
        ("read2", [(1000, 1100)], "100M"),                        # 单页
        ("read3", [(5990, 6040), (13050, 13100)], "50M7010N50M"),  # M 段自身跨页 0->1, 再 ->2
    ]

    refp = os.path.join(outdir, "ref.fa")
    with open(refp, "w") as f:
        f.write(">%s\n" % CHR)
        for i in range(0, LEN, LINE):
            f.write(seq[i:i + LINE] + "\n")

    samp = os.path.join(outdir, "case.sam")
    expf = os.path.join(outdir, "expected_pos.list")
    expected = set()
    with open(samp, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n")
        f.write("@SQ\tSN:%s\tLN:%d\n" % (CHR, LEN))
        for name, wins, cigar in sorted(reads, key=lambda r: r[1][0][0]):
            rseq = "".join(seq[s:e] for s, e in wins)
            pos = wins[0][0] + 1  # 1-based
            f.write(
                "%s\t0\t%s\t%d\t60\t%s\t*\t0\t0\t%s\t%s\t"
                "NM:i:0\tAS:i:%d\tNH:i:1\tMD:Z:%d\tYZ:A:+\tXN:i:0\n"
                % (name, CHR, pos, cigar, rseq, "I" * len(rseq),
                   len(rseq), len(rseq)))
            for s, e in wins:
                for i in range(s, e):
                    if seq[i] == 'C':
                        expected.add(i + 1)
    with open(expf, "w") as f:
        for p in sorted(expected):
            f.write("%d\n" % p)

    print("ref   : %s (%d bp)" % (refp, LEN))
    print("sam   : %s (%d reads)" % (samp, len(reads)))
    print("oracle: %s (%d positions)" % (expf, len(expected)))


if __name__ == "__main__":
    sys.exit(main())
