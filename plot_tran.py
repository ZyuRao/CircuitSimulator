#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
读取瞬态结果 CSV（tran_out.csv），画时域波形。

用法：
    python plot_tran.py tran_out.csv 'V(118)'
    python plot_tran.py tran_out.csv 'V(101)' 'V(118)'
    python plot_tran.py tran_out.csv    # 自动画所有电压列
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_tran.py <csv_file> [col1] [col2] ... [--out=png]")
        sys.exit(1)

    csv_file = sys.argv[1]
    out_file = None

    cols = []
    args = sys.argv[2:]
    i = 0
    while i < len(args):
        a = args[i]
        if a.startswith("--out="):
            out_file = a.split("=", 1)[1]
        elif a == "--out":
            if i + 1 < len(args):
                out_file = args[i + 1]
                i += 1
            else:
                print("--out requires a path")
                sys.exit(1)
        else:
            cols.append(a)
        i += 1

    if not os.path.exists(csv_file):
        print(f"File not found: {csv_file}")
        sys.exit(1)

    # 读取 CSV
    df = pd.read_csv(csv_file)

    if "time" not in df.columns:
        print("CSV 文件中没有 'time' 列，确认一下输出格式。")
        print("当前列名：", list(df.columns))
        sys.exit(1)

    t = df["time"]

    # 要画哪些列？
    if not cols:
        # 没指定时：自动选所有电压列（列名形如 V(XXX)）
        cols = [c for c in df.columns if c.startswith("V(")]
        if not cols:
            print("未找到任何以 'V(' 开头的电压列，请手动指定要画的列名。")
            print("当前列名：", list(df.columns))
            sys.exit(1)
        print("自动选择的电压列：", cols)

    # 检查列是否存在
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print("以下列名在 CSV 中不存在：", missing)
        print("当前列名：", list(df.columns))
        sys.exit(1)

    # 画图
    plt.figure()
    for c in cols:
        plt.plot(t, df[c], label=c)

    plt.xlabel("Time (s)")
    plt.ylabel("Value")
    plt.title(os.path.basename(csv_file))
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    if out_file:
        os.makedirs(os.path.dirname(out_file), exist_ok=True)
        plt.savefig(out_file)
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    main()
