#!/usr/bin/env python
# -*- coding=utf-8 -*-
"""
TMB Cohort Summary v1
======================
Aggregate per-sample TMB result files into a cohort-level summary.

Input: a directory containing per-sample *.tmb.tsv files.
输出：
  - cohort_tmb_summary.tsv   : 所有样本横表（一行一样本）
  - cohort_tmb_summary.txt   : 队列层统计（中位数、IQR、TMB-High 比例等）
  - cohort_tmb_distribution.png  : TMB 分布图（如果装了 matplotlib）

用法：
  python summarize_tmb.py -d /path/to/tmb_dir -o /path/to/output_dir
  python summarize_tmb.py -d /path/to/tmb_dir -o /path/to/output_dir --label COHORT001
  python summarize_tmb.py -d /path/to/tmb_dir -o /path/to/output_dir --high-cutoff 10
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


def get_opt():
    p = argparse.ArgumentParser(
        description="Aggregate per-sample TMB results into cohort summary",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-d", "--tmb-dir", required=True, metavar="DIR",
                   help="Directory containing per-sample *.tmb.tsv files")
    p.add_argument("-o", "--out-dir", required=True, metavar="DIR",
                   help="Output directory for cohort summary")
    p.add_argument("--label", default="cohort", metavar="STR",
                   help="Cohort label (used in filenames and titles)")
    p.add_argument("--high-cutoff", type=float, default=10.0, metavar="FLOAT",
                   help="TMB-High cutoff [default: 10.0]")
    p.add_argument("--pattern", default="*.tmb.tsv", metavar="GLOB",
                   help="Glob pattern for per-sample files [default: *.tmb.tsv]")
    p.add_argument("--no-plot", action="store_true",
                   help="Skip distribution plot")
    return p.parse_args()


def main():
    opt = get_opt()

    tmb_dir = Path(opt.tmb_dir)
    out_dir = Path(opt.out_dir)
    if not tmb_dir.is_dir():
        sys.exit(f"Error: TMB directory not found: {tmb_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # ===== 收集所有文件 =====
    files = sorted(tmb_dir.glob(opt.pattern))
    if not files:
        sys.exit(f"Error: no files matching {opt.pattern} in {tmb_dir}")

    print(f"\n{'=' * 70}")
    print(f"TMB Cohort Summary: {opt.label}")
    print(f"{'=' * 70}")
    print(f"Input dir:    {tmb_dir}")
    print(f"Pattern:      {opt.pattern}")
    print(f"Files found:  {len(files)}")
    print(f"Output dir:   {out_dir}")
    print(f"TMB-High:     ≥ {opt.high_cutoff} muts/Mb")
    print(f"{'=' * 70}\n")

    # ===== 拼接 =====
    rows = []
    failed = []
    for fp in files:
        try:
            df = pd.read_csv(fp, sep="\t")
            if len(df) == 0:
                failed.append((fp.name, "empty file"))
                continue
            # 单样本文件预期 1 行；多样本也兼容
            for _, r in df.iterrows():
                rows.append(r)
        except Exception as e:
            failed.append((fp.name, str(e)))

    if not rows:
        sys.exit("Error: no usable rows in any TMB file")

    cohort = pd.DataFrame(rows).reset_index(drop=True)
    print(f"Successfully loaded {len(cohort)} samples "
          f"({len(failed)} failed)\n")

    # ===== 写横表 =====
    table_path = out_dir / f"{opt.label}_tmb_summary.tsv"
    cohort.to_csv(table_path, sep="\t", index=False)
    print(f"✓ Sample-level table:  {table_path}")

    # ===== 队列统计 =====
    if "TMB_Nonsynonymous" not in cohort.columns:
        sys.exit("Error: column TMB_Nonsynonymous missing in input files")

    tmb = pd.to_numeric(cohort["TMB_Nonsynonymous"], errors="coerce").dropna()
    n = len(tmb)
    n_high = (tmb >= opt.high_cutoff).sum()

    stats_lines = [
        f"# TMB Cohort Summary: {opt.label}",
        f"# Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"# Source: {tmb_dir}",
        "",
        f"Total samples loaded:      {len(cohort)}",
        f"Samples with valid TMB:    {n}",
        f"Failed/empty samples:      {len(failed)}",
        "",
        "==== TMB_Nonsynonymous (FoCR-aligned) ====",
        f"  Mean:        {tmb.mean():.3f} muts/Mb",
        f"  Median:      {tmb.median():.3f} muts/Mb",
        f"  Std dev:     {tmb.std():.3f}",
        f"  Min:         {tmb.min():.3f}",
        f"  Q1 (25%):    {tmb.quantile(0.25):.3f}",
        f"  Q3 (75%):    {tmb.quantile(0.75):.3f}",
        f"  Max:         {tmb.max():.3f}",
        f"  IQR:         {tmb.quantile(0.75) - tmb.quantile(0.25):.3f}",
        "",
        f"==== TMB-High proportion (cutoff = {opt.high_cutoff}) ====",
        f"  TMB-High:    {n_high}/{n} ({n_high/n*100:.2f}%)",
        f"  TMB-Low:     {n-n_high}/{n} ({(n-n_high)/n*100:.2f}%)",
        "",
        "==== Distribution buckets ====",
    ]
    bins = [
        ("< 1",        (tmb < 1).sum()),
        ("1 - 3",      ((tmb >= 1) & (tmb < 3)).sum()),
        ("3 - 6",      ((tmb >= 3) & (tmb < 6)).sum()),
        ("6 - 10",     ((tmb >= 6) & (tmb < 10)).sum()),
        ("10 - 20",    ((tmb >= 10) & (tmb < 20)).sum()),
        ("20 - 50",    ((tmb >= 20) & (tmb < 50)).sum()),
        ("≥ 50",       (tmb >= 50).sum()),
    ]
    for lab, cnt in bins:
        stats_lines.append(f"  TMB {lab:10s} muts/Mb: {cnt:4d} ({cnt/n*100:5.1f}%)")

    if failed:
        stats_lines += ["", "==== Failed files ===="]
        for name, err in failed[:20]:
            stats_lines.append(f"  {name}: {err}")
        if len(failed) > 20:
            stats_lines.append(f"  ... and {len(failed)-20} more")

    stats_text = "\n".join(stats_lines)

    stats_path = out_dir / f"{opt.label}_tmb_summary.txt"
    with open(stats_path, "w") as f:
        f.write(stats_text + "\n")
    print(f"✓ Cohort statistics:   {stats_path}")
    print()
    print(stats_text)

    # ===== 分布图 =====
    if not opt.no_plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # 左：histogram (linear scale, capped at 95th percentile)
            ax = axes[0]
            cap = tmb.quantile(0.95) * 1.5 if tmb.quantile(0.95) > 0 else tmb.max()
            ax.hist(tmb[tmb <= cap], bins=40, edgecolor="black", alpha=0.7,
                    color="steelblue")
            ax.axvline(opt.high_cutoff, color="red", linestyle="--", linewidth=2,
                       label=f"TMB-High cutoff ({opt.high_cutoff})")
            ax.axvline(tmb.median(), color="orange", linestyle="-", linewidth=2,
                       label=f"Median ({tmb.median():.2f})")
            ax.set_xlabel("TMB (muts/Mb)")
            ax.set_ylabel("Number of samples")
            ax.set_title(f"{opt.label} — TMB distribution (n={n})")
            ax.legend()
            ax.grid(alpha=0.3)

            # 右：log-scale histogram (full range)
            ax = axes[1]
            tmb_for_log = tmb[tmb > 0]
            if len(tmb_for_log) > 0:
                bins_log = np.logspace(np.log10(max(tmb_for_log.min(), 0.01)),
                                       np.log10(tmb_for_log.max()), 40)
                ax.hist(tmb_for_log, bins=bins_log, edgecolor="black", alpha=0.7,
                        color="darkseagreen")
                ax.set_xscale("log")
                ax.axvline(opt.high_cutoff, color="red", linestyle="--", linewidth=2,
                           label=f"TMB-High ({opt.high_cutoff})")
                ax.set_xlabel("TMB (muts/Mb, log scale)")
                ax.set_ylabel("Number of samples")
                ax.set_title(f"{opt.label} — TMB distribution (log)")
                ax.legend()
                ax.grid(alpha=0.3, which="both")

            plt.tight_layout()
            plot_path = out_dir / f"{opt.label}_tmb_distribution.png"
            plt.savefig(plot_path, dpi=120)
            plt.close()
            print(f"✓ Distribution plot:   {plot_path}")
        except ImportError:
            print("⚠️  matplotlib not installed, skipping plot")
        except Exception as e:
            print(f"⚠️  plot generation failed: {e}")

    print(f"\n{'=' * 70}\n")


if __name__ == "__main__":
    main()
