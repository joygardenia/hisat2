#!/usr/bin/env bash
# 跨页读段回归测试。
#
# 用法:
#   make hisat_3n_table
#   bash tests/cross_page/run.sh [binary]
#
# 断言:
#   1. 程序退出码为 0(不崩溃)
#   2. 输出中 strand='+' 的位点集合 == expected_pos.list(oracle)
#
# 修复前:退出码 134(越界 -> std::system_error),失败。
# 修复后:退出码 0 且位点集合一致,通过。
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="${1:-$HERE/../../hisat-3n-table}"
WORK="${WORK:-$HERE}"

if [ ! -x "$BIN" ]; then
  echo "FAIL: 找不到可执行文件 $BIN (先 make hisat_3n_table)" >&2
  exit 2
fi

python3 "$HERE/gen_fixture.py" "$WORK" || exit 2
samtools faidx "$WORK/ref.fa" || exit 2

OUT="$WORK/out.tsv"
"$BIN" --alignments "$WORK/case.sam" --ref "$WORK/ref.fa" \
       --output-name "$OUT" --base-change C,T
rc=$?
if [ $rc -ne 0 ]; then
  echo "FAIL: 退出码 $rc (期望 0) —— 跨页越界崩溃?" >&2
  exit 1
fi

python3 - "$OUT" "$WORK/expected_pos.list" <<'PY'
import sys

out_path, exp_path = sys.argv[1], sys.argv[2]
got = set()
with open(out_path) as f:
    for line in f:
        if line.startswith("ref\t"):
            continue
        col = line.rstrip("\n").split("\t")
        if len(col) >= 3 and col[2] == "+":
            got.add(int(col[1]))

want = set(int(x) for x in open(exp_path) if x.strip())

if got != want:
    print("FAIL: 位点集合不一致 got=%d want=%d" % (len(got), len(want)))
    print("  多出:", sorted(got - want)[:10])
    print("  缺少:", sorted(want - got)[:10])
    sys.exit(1)

print("PASS: %d 个位点与 oracle 一致" % len(got))
PY
