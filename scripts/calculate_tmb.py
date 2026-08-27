#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
TMB calculation from an annotated MAF.

The script filters variants using configurable depth/VAF/population-frequency
criteria and reports nonsynonymous and nonsilent TMB. A BED file can be used
to calculate the callable territory; otherwise a configurable exome-size
denominator is used.

Examples:
  # Fixed exome-size denominator
  python calculate_tmb.py -i sample.raw.maf -o sample.tmb.tsv

  # 用实际 BED 做分母（推荐用 GENCODE CDS 或 callable BED）
  python calculate_tmb.py -i sample.raw.maf -o sample.tmb.tsv \\
      --bed /path/to/cds.bed

  # 保留中间过滤 MAF（debug 用）
  python calculate_tmb.py -i sample.raw.maf -o sample.tmb.tsv \\
      --save-filtered sample.tmb_filtered.maf -s sample.tmb_filter.stat
"""

import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from collections import OrderedDict

# ====================================================================
# 变异类型定义
# ====================================================================

NONSYNONYMOUS_VARTYPE = [
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Nonstop_Mutation",
    "Translation_Start_Site",     # 旧脚本缺
    "Splice_Site",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
]

SILENT_AND_NONCODING = [
    "Silent",
    "Intron",
    "5'UTR", "3'UTR",
    "5'Flank", "3'Flank",
    "IGR",
    "RNA",
    "Splice_Region",
    "Targeted_Region",
]

# Default WES denominator; override with --exome-size-mb or --bed
DEFAULT_EXOME_SIZE_MB = 38.0
TMB_HIGH_CUTOFF = 10.0


# ====================================================================
# 命令行
# ====================================================================

def get_opt():
    p = argparse.ArgumentParser(
        description="TMB Pipeline v4: filter + calculate TMB in one step",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # 输入输出
    p.add_argument("-i", required=True, metavar="FILE",
                   help="Input raw MAF (from vcf2maf)")
    p.add_argument("-o", required=True, metavar="FILE",
                   help="Output TMB result table (tsv)")
    p.add_argument("--save-filtered", default=None, metavar="FILE",
                   help="Optional: save the intermediate filtered MAF "
                        "(useful for debugging, default: not saved)")
    p.add_argument("-s", "--stat", default=None, metavar="FILE",
                   help="Optional: save filter step counts (default: not saved)")

    # 过滤参数（对齐国际标准的默认值）
    p.add_argument("-t_AF", type=float, default=0.05, metavar="FLOAT",
                   help="Min tumor VAF [default: 0.05, AMP 2024]")
    p.add_argument("-t_depth", type=int, default=25, metavar="INT",
                   help="Min tumor depth [default: 25x, AMP 2024]")
    p.add_argument("-n_depth", type=int, default=10, metavar="INT",
                   help="Min normal depth [default: 10x, MC3]")
    p.add_argument("-t_alt_count", type=int, default=3, metavar="INT",
                   help="Min tumor alt reads [default: 3]")
    p.add_argument("-n_alt_count", type=int, default=1, metavar="INT",
                   help="Max normal alt reads [default: 1, MC3/Yang 2024]")
    p.add_argument("-gnomAD_AF", type=float, default=0.001, metavar="FLOAT",
                   help="Max gnomAD population AF [default: 0.001]")

    # 分母选择
    p.add_argument("--exome-size-mb", type=float, default=DEFAULT_EXOME_SIZE_MB,
                   metavar="FLOAT",
                   help=f"Exome size in Mb [default: {DEFAULT_EXOME_SIZE_MB}, FoCR consensus]")
    p.add_argument("--bed", default=None, metavar="FILE",
                   help="Optional BED file: use sum of intervals as denominator "
                        "(more accurate than fixed 38 Mb)")

    # TMB-High cutoff
    p.add_argument("--high-cutoff", type=float, default=TMB_HIGH_CUTOFF,
                   metavar="FLOAT",
                   help=f"TMB-High cutoff [default: {TMB_HIGH_CUTOFF}, FDA]")

    return p.parse_args()


# ====================================================================
# 工具函数
# ====================================================================

def find_gnomad_column(maf):
    """自动识别 gnomAD AF 列（vcf2maf 不同版本列名可能不同）"""
    candidates = [c for c in maf.columns if "gnomad" in c.lower() and "af" in c.lower()]
    if not candidates:
        return None
    for kw in ["gnomAD_AF", "gnomad_af", "gnomADg_AF"]:
        for c in candidates:
            if c.lower() == kw.lower():
                return c
    return candidates[0]


def recalculate_af(maf):
    """重算 t_AF 和 n_AF；n_AF=0 是合法值（最干净的 somatic）"""
    maf = maf.copy(deep=True)
    for col in ["n_depth", "n_alt_count", "t_depth", "t_alt_count"]:
        if col not in maf.columns:
            raise KeyError(f"Missing required column: {col}")
        maf[col] = pd.to_numeric(maf[col], errors="coerce")

    with np.errstate(divide="ignore", invalid="ignore"):
        maf["t_AF"] = np.where(maf["t_depth"] > 0,
                               maf["t_alt_count"] / maf["t_depth"],
                               np.nan)
        maf["n_AF"] = np.where(maf["n_depth"] > 0,
                               maf["n_alt_count"] / maf["n_depth"],
                               np.nan)
    return maf


def get_denominator_mb(opt):
    """返回分母 (Mb) 和说明"""
    if opt.bed:
        bed_path = Path(opt.bed)
        if not bed_path.exists():
            sys.exit(f"Error: BED file not found: {opt.bed}")
        total_bp = 0
        with open(bed_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    try:
                        total_bp += int(fields[2]) - int(fields[1])
                    except ValueError:
                        continue
        mb = total_bp / 1e6
        return mb, f"BED file ({bed_path.name}, {mb:.2f} Mb)"
    return opt.exome_size_mb, f"Fixed {opt.exome_size_mb} Mb (FoCR WES consensus)"


# ====================================================================
# 主流程
# ====================================================================

def filter_maf(maf, opt, count):
    """对 raw MAF 应用 TMB-specific 过滤；返回过滤后的 MAF 和 counts"""

    # Step 1: FILTER=PASS
    if "FILTER" in maf.columns:
        maf = maf[maf["FILTER"] == "PASS"]
    else:
        print("  Warning: FILTER column missing, skipping FILTER=PASS step")
    count["FILTER_PASS"] = len(maf)
    print(f"  Step 1 [FILTER=PASS]:        {len(maf):,}")
    if maf.empty:
        return maf, count

    # Step 2: 重算 AF（修复 n_AF=0 安全处理）
    maf = recalculate_af(maf)

    # Step 3: 变异类型黑名单（TMB 用 nonsilent 定义）
    if "Variant_Classification" in maf.columns:
        # 先按黑名单过滤（去掉 Silent + 非编码）
        maf = maf[~maf["Variant_Classification"].isin(SILENT_AND_NONCODING)]
        count["VarType_nonsilent"] = len(maf)
        print(f"  Step 2 [non-silent (blacklist)]: {len(maf):,}")
    else:
        print("  Warning: Variant_Classification column missing, skipping")

    if maf.empty:
        return maf, count

    # Step 4-7: 深度 / alt count 过滤
    for label, expr, key in [
        (f"t_depth ≥ {opt.t_depth}",      maf["t_depth"]      >= opt.t_depth,      f"t_depth>={opt.t_depth}"),
        (f"n_depth ≥ {opt.n_depth}",      maf["n_depth"]      >= opt.n_depth,      f"n_depth>={opt.n_depth}"),
        (f"t_alt   ≥ {opt.t_alt_count}",  maf["t_alt_count"]  >= opt.t_alt_count,  f"t_alt>={opt.t_alt_count}"),
        (f"n_alt   ≤ {opt.n_alt_count}",  maf["n_alt_count"]  <= opt.n_alt_count,  f"n_alt<={opt.n_alt_count}"),
        (f"t_AF    ≥ {opt.t_AF}",         maf["t_AF"]         >= opt.t_AF,         f"t_AF>={opt.t_AF}"),
    ]:
        maf = maf[expr]
        count[key] = len(maf)
        print(f"  Step [{label}]: {len(maf):,}")
        if maf.empty:
            return maf, count

    # Step 8: gnomAD 群体频率
    gnomad_col = find_gnomad_column(maf)
    if gnomad_col is not None:
        af = pd.to_numeric(maf[gnomad_col], errors="coerce")
        keep = af.isna() | (af < opt.gnomAD_AF)
        maf = maf[keep]
        count[f"gnomAD<{opt.gnomAD_AF}"] = len(maf)
        print(f"  Step [gnomAD AF < {opt.gnomAD_AF}]: {len(maf):,}  (using col: {gnomad_col})")
    else:
        print("  Warning: no gnomAD column found, skipping population AF filter")

    return maf, count


def calc_tmb(maf, denom_mb, high_cutoff):
    """对（已过滤的）MAF 按样本计算 TMB"""

    if "Variant_Classification" in maf.columns:
        n_nonsyn = maf["Variant_Classification"].isin(NONSYNONYMOUS_VARTYPE).sum()
        n_nonsilent = (~maf["Variant_Classification"].isin(SILENT_AND_NONCODING)).sum()
    else:
        n_nonsyn = len(maf)
        n_nonsilent = len(maf)

    tmb_nonsyn = round(n_nonsyn / denom_mb, 4) if denom_mb > 0 else 0
    tmb_nonsilent = round(n_nonsilent / denom_mb, 4) if denom_mb > 0 else 0

    avg_t_depth = "NA"
    if "t_depth" in maf.columns and len(maf) > 0:
        d = pd.to_numeric(maf["t_depth"], errors="coerce")
        if d.notna().any():
            avg_t_depth = int(d.mean())

    return {
        "Mutation_Count_Total": len(maf),
        "Mutation_Count_Nonsynonymous_9types": int(n_nonsyn),
        "Mutation_Count_Nonsilent": int(n_nonsilent),
        "TMB_Nonsynonymous": tmb_nonsyn,
        "TMB_Nonsilent": tmb_nonsilent,
        "TMB_Category": "High" if tmb_nonsyn >= high_cutoff else "Low",
        "Average_Tumor_Depth": avg_t_depth,
    }


def write_empty_tmb(opt, denom_mb, sample_id):
    """空输入也写一份带表头的结果，避免下游 pipeline 失败"""
    df = pd.DataFrame([{
        "Sample_ID": sample_id,
        "Mutation_Count_Total": 0,
        "Mutation_Count_Nonsynonymous_9types": 0,
        "Mutation_Count_Nonsilent": 0,
        "Exome_Size_Mb": round(denom_mb, 2),
        "TMB_Nonsynonymous": 0.0,
        "TMB_Nonsilent": 0.0,
        "TMB_Category": "Low",
        "Average_Tumor_Depth": "NA",
    }])
    df.to_csv(opt.o, sep="\t", index=False)


def main():
    opt = get_opt()

    # 检查输入
    in_path = Path(opt.i)
    if not in_path.exists():
        sys.exit(f"Error: Input MAF not found: {opt.i}")
    if in_path.suffix != ".maf":
        print(f"Warning: Input doesn't have .maf extension: {opt.i}")

    print(f"\n{'=' * 72}")
    print(f"TMB Pipeline v4 — Filter + Calculate (FoCR/AMP/MC3-aligned)")
    print(f"{'=' * 72}")
    print(f"Input MAF:    {opt.i}")
    print(f"Output TMB:   {opt.o}")

    denom_mb, denom_note = get_denominator_mb(opt)
    print(f"Denominator:  {denom_note}")
    print(f"TMB-High:     ≥ {opt.high_cutoff} muts/Mb (FDA pembrolizumab)")
    print(f"Filter params:")
    print(f"  t_AF≥{opt.t_AF}  t_depth≥{opt.t_depth}  n_depth≥{opt.n_depth}  "
          f"t_alt≥{opt.t_alt_count}  n_alt≤{opt.n_alt_count}  gnomAD<{opt.gnomAD_AF}")
    print(f"{'=' * 72}\n")

    # 读 MAF
    try:
        maf = pd.read_csv(opt.i, sep="\t", comment="#", low_memory=False)
    except Exception as e:
        sys.exit(f"Error reading MAF: {e}")

    count = OrderedDict()
    count["Total"] = len(maf)
    print(f"Loaded {len(maf):,} variants from MAF\n")

    # 推断 sample_id（优先用 MAF 内的 Tumor_Sample_Barcode；否则用文件名）
    if "Tumor_Sample_Barcode" in maf.columns and len(maf) > 0:
        unique_samples = maf["Tumor_Sample_Barcode"].dropna().unique()
    else:
        unique_samples = []
    fallback_id = in_path.stem.replace(".vep", "").replace("_filtered", "")

    # 空 MAF
    if len(maf) == 0:
        sample_id = fallback_id
        print(f"Empty MAF — writing zero TMB for {sample_id}")
        write_empty_tmb(opt, denom_mb, sample_id)
        return

    # 过滤
    print("Filtering MAF for TMB:")
    filtered, count = filter_maf(maf, opt, count)
    print(f"\nFinal filtered variants: {len(filtered):,} "
          f"({len(filtered) / count['Total'] * 100:.1f}% of input)\n")

    # 可选：保存中间 MAF
    if opt.save_filtered:
        if filtered.empty:
            # 写表头
            with open(opt.save_filtered, "w") as fo:
                with open(opt.i) as fi:
                    for line in fi:
                        if line.startswith("Hugo_Symbol"):
                            fo.write(line)
                            break
        else:
            filtered.to_csv(opt.save_filtered, sep="\t", index=False)
        print(f"✓ Intermediate filtered MAF saved: {opt.save_filtered}")

    # 可选：保存过滤步骤统计
    if opt.stat:
        df_count = pd.DataFrame.from_dict(count, orient="index", columns=["Variants"])
        df_count.to_csv(opt.stat, sep="\t", header=False)
        print(f"✓ Filter stats saved: {opt.stat}")

    # 计算 TMB（按样本分组；通常单样本 MAF 只有一个）
    if filtered.empty:
        sample_id = unique_samples[0] if len(unique_samples) > 0 else fallback_id
        print(f"All variants filtered out — TMB = 0 for {sample_id}")
        write_empty_tmb(opt, denom_mb, sample_id)
        return

    if "Tumor_Sample_Barcode" in filtered.columns:
        results = []
        for sample_id, sub in filtered.groupby("Tumor_Sample_Barcode"):
            r = calc_tmb(sub, denom_mb, opt.high_cutoff)
            r["Sample_ID"] = sample_id
            r["Exome_Size_Mb"] = round(denom_mb, 2)
            results.append(r)

            print(f"Sample {sample_id}:")
            print(f"  Variants kept:           {r['Mutation_Count_Total']:,}")
            print(f"  Nonsynonymous (9 types): {r['Mutation_Count_Nonsynonymous_9types']:,}")
            print(f"  Nonsilent:               {r['Mutation_Count_Nonsilent']:,}")
            print(f"  TMB_Nonsynonymous:       {r['TMB_Nonsynonymous']:.2f} muts/Mb  [{r['TMB_Category']}]")
            print(f"  TMB_Nonsilent:           {r['TMB_Nonsilent']:.2f} muts/Mb")
            print(f"  Mean tumor depth:        {r['Average_Tumor_Depth']}")
    else:
        sample_id = fallback_id
        r = calc_tmb(filtered, denom_mb, opt.high_cutoff)
        r["Sample_ID"] = sample_id
        r["Exome_Size_Mb"] = round(denom_mb, 2)
        results = [r]

        print(f"Sample {sample_id}:")
        print(f"  TMB_Nonsynonymous:       {r['TMB_Nonsynonymous']:.2f} muts/Mb  [{r['TMB_Category']}]")
        print(f"  TMB_Nonsilent:           {r['TMB_Nonsilent']:.2f} muts/Mb")

    # 输出
    cols = ["Sample_ID",
            "Mutation_Count_Total",
            "Mutation_Count_Nonsynonymous_9types",
            "Mutation_Count_Nonsilent",
            "Exome_Size_Mb",
            "TMB_Nonsynonymous",
            "TMB_Nonsilent",
            "TMB_Category",
            "Average_Tumor_Depth"]
    df = pd.DataFrame(results)[cols]
    df.to_csv(opt.o, sep="\t", index=False)

    print(f"\n✓ TMB result saved: {opt.o}")
    print(f"{'=' * 72}\n")


if __name__ == "__main__":
    main()
