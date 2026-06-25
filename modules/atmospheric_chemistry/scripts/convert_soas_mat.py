#!/usr/bin/env python3
"""
convert_soas_mat.py — 将 F0AM SOAS 观测 .mat 文件转换为 MOOSE 可读的 CSV 格式。

用法:
    python3 convert_soas_mat.py <mat_file> [--output <csv_file>]

参数:
    mat_file        SOAS 观测 .mat 文件路径（必需）
    --output, -o    输出 CSV 文件路径（默认: 与输入同目录下的 SOAS_DielCycle.csv）
    --hour, -H      仅输出指定小时的初始条件（默认: 全部 24 小时）

示例:
    # 转换全部 24 小时数据
    python3 convert_soas_mat.py Obs_SOAS_CampaignAvg_60min.mat

    # 指定输出路径
    python3 convert_soas_mat.py Obs_SOAS_CampaignAvg_60min.mat -o ../database/SOAS_DielCycle.csv

    # 仅提取 hour 0 初始条件
    python3 convert_soas_mat.py Obs_SOAS_CampaignAvg_60min.mat --hour 0

输出 CSV 列:
    Time_h, T(K), P(mbar), RH(%), SZA(deg), M(molec/cm3), BLheight(m),
    NO, NO2, O3, OH, CO, H2O2, PAN, C2H5NO3, IC3H7NO3,
    C5H8, APINENE, BPINENE, LIMONENE, C2H4, C2H6, C3H8, ...
    H2 (常数 550 ppb), CH4 (常数 1770 ppb)

依赖: scipy, numpy
"""

import scipy.io as sio
import numpy as np
import argparse
import os
import sys

# SOAS 结构体中包含的观测物种字段（单位: ppb）
SPECIES_FIELDS = [
    'NO', 'NO2', 'O3', 'OH', 'CO', 'H2O2', 'PAN', 'C2H5NO3', 'IC3H7NO3',
    'C5H8', 'APINENE', 'BPINENE', 'LIMONENE',
    'C2H4', 'C2H6', 'C3H8', 'IC4H10', 'IC5H12', 'NC5H12', 'NC6H14', 'NC10H22',
    'BENZENE', 'TOLUENE', 'EBENZ', 'TM124B', 'TM135B', 'MXYL', 'OXYL', 'PXYL', 'BENZAL',
    'CH3CHO', 'C2H5CHO', 'C3H7CHO', 'HOCH2CHO', 'GLYOX', 'CH3OH', 'C2H5OH',
    'ACETOL', 'BIACET', 'MACR', 'MVK', 'HCHO', 'CH3COCH3',
    'C2H2', 'C3H6', 'NC4H10', 'DMS', 'HNO3', 'HO2', 'IEPOX', 'ISOPOOH', 'MEK', 'MPAN',
]

# 环境变量字段（原始单位）
ENV_FIELDS = ['T', 'P', 'RH', 'SZA', 'M', 'BLheight']

# 常数物种（不在 SOAS 观测中，在 F0AM 脚本中手动设定）
CONSTANT_SPECIES = {'H2': 550.0, 'CH4': 1770.0}  # ppb


def get_field_1d(soas, name):
    """从 SOAS 结构体中提取一维 numpy 数组。"""
    val = soas[name]
    if val.ndim == 0:
        val = val[()]
    return np.asarray(val, dtype=float).ravel()


def main():
    parser = argparse.ArgumentParser(
        description='将 F0AM SOAS 观测 .mat 文件转换为 MOOSE CSV 格式',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s Obs_SOAS_CampaignAvg_60min.mat
  %(prog)s Obs_SOAS_CampaignAvg_60min.mat -o SOAS_DielCycle.csv
  %(prog)s Obs_SOAS_CampaignAvg_60min.mat --hour 0
        """
    )
    parser.add_argument('mat_file', help='SOAS 观测 .mat 文件路径')
    parser.add_argument('--output', '-o', default=None,
                        help='输出 CSV 文件路径（默认: 与输入同目录下的 SOAS_DielCycle.csv）')
    parser.add_argument('--hour', '-H', type=int, default=None,
                        help='仅输出指定小时的初始条件')
    args = parser.parse_args()

    # 检查输入文件
    if not os.path.isfile(args.mat_file):
        print(f"错误: 找不到文件 '{args.mat_file}'", file=sys.stderr)
        sys.exit(1)

    # 默认输出路径
    out_path = args.output
    if out_path is None:
        out_dir = os.path.dirname(os.path.abspath(args.mat_file)) or '.'
        out_path = os.path.join(out_dir, 'SOAS_DielCycle.csv')

    # 读取 .mat 文件
    mat = sio.loadmat(args.mat_file, squeeze_me=True)
    if 'SOAS' not in mat:
        print("错误: .mat 文件中未找到 'SOAS' 结构体", file=sys.stderr)
        sys.exit(1)
    soas = mat['SOAS']

    time_arr = get_field_1d(soas, 'Time')
    n = len(time_arr)
    print(f"SOAS: {n} 个时间点, {len(soas.dtype.names)} 个字段")

    # 确定输出行范围
    if args.hour is not None:
        hour_indices = [args.hour]
    else:
        hour_indices = range(n)

    # 收集列名
    all_cols = ['Time_h'] + ENV_FIELDS + SPECIES_FIELDS + list(CONSTANT_SPECIES.keys())

    # 构建数据行
    rows = []
    for i in hour_indices:
        row = {'Time_h': int(time_arr[i])}
        for f in ENV_FIELDS:
            vals = get_field_1d(soas, f)
            row[f] = float(vals[i]) if i < len(vals) else 0.0
        for f in SPECIES_FIELDS:
            if f in soas.dtype.names:
                vals = get_field_1d(soas, f)
                row[f] = float(vals[i]) if i < len(vals) else 0.0
            else:
                row[f] = 0.0
        for name, val in CONSTANT_SPECIES.items():
            row[name] = val
        rows.append(row)

    # 确保输出目录存在
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)

    # 写 CSV
    with open(out_path, 'w') as f:
        f.write(','.join(all_cols) + '\n')
        for row in rows:
            f.write(','.join(str(row[c]) for c in all_cols) + '\n')

    print(f"已写入 {len(rows)} 行 × {len(all_cols)} 列 → {out_path}")

    # 打印初始条件摘要
    if rows:
        print(f"\n初始条件 (hour {rows[0]['Time_h']}):")
        for f in ENV_FIELDS:
            print(f"  {f}: {rows[0][f]:.6g}")
        print(f"  常数 H2: {CONSTANT_SPECIES['H2']} ppb")
        print(f"  常数 CH4: {CONSTANT_SPECIES['CH4']} ppb")
        for f in SPECIES_FIELDS[:8]:
            if f in rows[0]:
                print(f"  {f}: {rows[0][f]:.6g} ppb")
        print("  ...")


if __name__ == '__main__':
    main()
