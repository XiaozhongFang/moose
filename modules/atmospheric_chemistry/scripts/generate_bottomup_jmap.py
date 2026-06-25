#!/usr/bin/env python3
"""
generate_bottomup_jmap.py — 预计算 BottomUp 光解的截面/量子产率 CSV 并生成反应映射文件。

读取 F0AM 的 J_BottomUp.m 反应列表，解析每个光解反应的 CS/QY 数据源，
在指定温度下预计算截面和量子产率，输出统一的 2 列 CSV 和反应映射文件。

用法:
    python3 generate_bottomup_jmap.py <photolysis_dir> <output_dir> [T]

参数:
    photolysis_dir   F0AM Chem/Photolysis 目录路径（含 CrossSections/ 和 QuantumYields/）
    output_dir       输出目录（将在此创建 CrossSections_precomputed/ 等子目录）
    T                预计算温度 (K)，默认 298.0

示例:
    python3 generate_bottomup_jmap.py \\
        /path/to/F0AM/Chem/Photolysis \\
        /path/to/moose/database/photolysis/bottomup \\
        298.0

输出文件:
    output_dir/CrossSections/       — 2 列 CSV (波长nm, 截面cm²)
    output_dir/QuantumYields/       — 2 列 CSV (波长nm, 量子产率)
    output_dir/bottomup_jmap.dat    — 反应映射（JNAME CS_FILE CS_TYPE QY_FILE QY_TYPE）

映射文件格式:
    # 注释行
    JNAME  CS_FILE  CS_TYPE  QY_FILE  QY_TYPE
    CS_TYPE: 1=2列CSV, 2=3列CSV(T内插), 3=TXT
    QY_TYPE: 0=标量, 1=2列CSV, 2=3列CSV, 3=TXT

依赖: Python 3.8+, numpy (无其他依赖)
"""

import os, sys, re, csv, math
from pathlib import Path

# ── J_BottomUp.m 的反应列表 ──
# 每个条目: (J名称, CS来源, QY来源)
# 来源格式: '@函数名' (MATLAB函数句柄), '文件名.csv' (CSV), '文件名.txt' (TXT), '标量' (QY常量)

J_BOTTOMUP_REACTIONS = [
    # === MCM 标准光解 J1-J61 ===
    ('J1',  '@Cross_Section_O3_JPL',           '@Quantum_Yield_O3_O1D_JPL'),
    ('J2',  '@Cross_Section_O3_JPL',           '@Quantum_Yield_O3_O3P_JPL'),
    ('J3',  '@Cross_Section_H2O2',             '1'),
    ('J4',  '@Cross_Section_NO2',              '@Quantum_Yield_NO2'),
    ('J5',  '@Cross_Section_NO3',              'Quantum_Yield_NO3_NO.csv'),
    ('J6',  '@Cross_Section_NO3',              'Quantum_Yield_NO3_NO2.csv'),
    ('J7',  'Cross_Section_HONO.csv',          '1'),
    ('J8',  '@Cross_Section_HNO3',             '1'),
    ('J11', '@Cross_Section_HCHO',             '@Quantum_Yield_HCHO_HCO'),
    ('J12', '@Cross_Section_HCHO',             '@Quantum_Yield_HCHO_H2'),
    ('J13', 'Cross_Section_CH3CHO.csv',        '@Quantum_Yield_CH3CHO_CH3'),
    ('J14', 'Cross_Section_C2H5CHO.csv',        'Quantum_Yield_C2H5CHO.csv'),
    ('J15', 'Cross_Section_C3H7CHO.csv',        '@Quantum_Yield_C3H7CHO_C3H7'),
    ('J16', 'Cross_Section_C3H7CHO.csv',        '@Quantum_Yield_C3H7CHO_C2H4'),
    ('J17', 'Cross_Section_IPRCHO.csv',         'Quantum_Yield_IPRCHO.csv'),
    ('J18', 'Cross_Section_MACR.csv',           '0.005'),
    ('J19', 'Cross_Section_MACR.csv',           '0.005'),
    ('J20', 'Cross_Section_MACR.csv',           '1'),
    ('J21', '@Cross_Section_CH3COCH3',         '@Quantum_Yield_CH3COCH3_CH3CO'),
    ('J22', 'Cross_Section_MEK.csv',            '0.34'),
    ('J23', 'Cross_Section_MVK.csv',           '@Quantum_Yield_MVK'),
    ('J24', 'Cross_Section_MVK.csv',           '@Quantum_Yield_MVK'),
    ('J31', 'Cross_Section_GLYOX.csv',         '@Quantum_Yield_GLYOX_H2'),
    ('J32', 'Cross_Section_GLYOX.csv',         '@Quantum_Yield_GLYOX_HCHO'),
    ('J33', 'Cross_Section_GLYOX.csv',         '@Quantum_Yield_GLYOX_HCO'),
    ('J34', 'Cross_Section_MGLYOX.csv',        '@Quantum_Yield_MGLYOX'),
    ('J35', 'Cross_Section_BIACET.csv',         '0.158'),
    ('J41', 'Cross_Section_CH3OOH.csv',         '1'),
    ('J51', '@Cross_Section_CH3NO3',           '1'),
    ('J52', '@Cross_Section_C2H5NO3',          '1'),
    ('J53', 'Cross_Section_NC3H7NO3.csv',       '1'),
    ('J54', '@Cross_Section_IC3H7NO3',         '1'),
    ('J55', 'Cross_Section_TC4H9NO3.csv',       '1'),
    ('J56', 'Cross_Section_NOA.csv',            '0.9'),
    ('J57', 'Cross_Section_NOA.csv',            '0.1'),
    # === 非 MCM 光解 Jn1-Jn56 ===
    ('Jn1',  'Cross_Section_CH3CHCHCHO.csv',    '0.030'),
    ('Jn2',  'Cross_Section_C6H5COH.csv',       '0.29'),
    ('Jn3',  'Cross_Section_C2H5COC2H5.csv',    '1'),
    ('Jn4',  'Cross_Section_CH3COOOH.csv',      '1'),
    ('Jn5',  'Cross_Section_CH3CHO.csv',        '@Quantum_Yield_CH3CHO_CH4'),
    ('Jn6',  'Cross_Section_CH3CHO.csv',        '@Quantum_Yield_CH3CHO_CH3CO'),
    ('Jn8',  '@Cross_Section_CH3COCH3',        '@Quantum_Yield_CH3COCH3_CO'),
    ('Jn9',  'Cross_Section_GLYCOALDEHYDE.csv', '1'),
    ('Jn10', 'Cross_Section_Hydroxyaceton.csv', '0.60'),
    ('Jn11', 'Cross_Section_Acrolein.csv',     '@Quantum_Yield_Acrolein'),
    ('Jn12', 'Cross_Section_3_methyl_2_nitrophenol.csv', '0.00015'),
    ('Jn13', 'Cross_Section_4_methyl_2_nitrophenol.csv', '0.0001'),
    ('Jn14', '@Cross_Section_PAN',             'Quantum_Yield_PAN_NO2.csv'),
    ('Jn15', '@Cross_Section_PAN',             'Quantum_Yield_PAN_NO3.csv'),
    ('Jn16', 'Cross_Section_CH3O2NO2.csv',      '0.95'),
    ('Jn17', 'Cross_Section_CH3O2NO2.csv',      '0.05'),
    ('Jn18', 'Cross_Section_CH3ONO.csv',        '0.76'),
    ('Jn19', '@Cross_Section_N2O5',            'Quantum_Yield_N2O5_NO3_NO2.csv'),
    ('Jn20', '@Cross_Section_N2O5',            'Quantum_Yield_N2O5_NO3_NO_O.csv'),
    ('Jn21', 'Cross_Section_HO2NO2.csv',        '0.59'),
    ('Jn22', 'Cross_Section_HO2NO2.csv',        '0.41'),
    ('Jn23', '@Cross_Section_ClNO2',           '1'),
    ('Jn24', '@Cross_Section_Br2',             '1'),
    ('Jn25', 'Cross_Section_BrO.csv',           '1'),
    ('Jn26', 'Cross_Section_HOBr.csv',          '1'),
    ('Jn27', 'Cross_Section_BrNO2.csv',         '1'),
    ('Jn28', '@Cross_Section_BrONO2',          '0.85'),
    ('Jn29', '@Cross_Section_BrONO2',          '0.15'),
    ('Jn30', '@Cross_Section_CHBr3',           '1'),
    ('Jn31', 'Cross_Section_BrCl.csv',          '1'),
    ('Jn32', '@Cross_Section_Cl2',             '1'),
    ('Jn33', '@Cross_Section_ClO_MB1999',      '1'),
    ('Jn34', '@Cross_Section_ClONO2',          'Quantum_Yield_ClONO2_Cl.csv'),
    ('Jn35', '@Cross_Section_ClONO2',          'Quantum_Yield_ClONO2_ClO.csv'),
    ('Jn36', 'Cross_Section_HOCl.csv',          '1'),
    ('Jn37', 'Cross_Section_OClO_Wahner(1987)_296K_245-475nm(0.22nm).txt', '1'),
    ('Jn38', 'Cross_Section_ClOOCl_JPL-2010(2011)_190-250K_200-420nm(rec).txt', '1'),
    ('Jn39', 'Cross_Section_ClOO_JPL-2010(2011)_191K_220-280nm(rec).txt', '1'),
    ('Jn40', 'Cross_Section_I2_JPL-2010(2011)_295K_185-700nm(rec).txt', 'QY_I2.txt'),
    ('Jn41', 'Cross_Section_HOI_JPL-2010(2011)_295-298K_280-480nm(rec).txt', '1'),
    ('Jn42', 'Cross_Section_IO_JPL-2010(2011)_298K_339-417nm(rec).txt', '0.91'),
    ('Jn43', 'Cross_Section_OIO_JPL-2010(2011)_295K_516-572nm(rec).txt', '1'),
    ('Jn44', 'Cross_Section_INO_JPL-2010(2011)_298K_223-460nm(rec).txt', '1'),
    ('Jn45', 'Cross_Section_INO2_JPL-2010(2011)_298K_210-380nm(rec).txt', '1'),
    ('Jn46', 'Cross_Section_IONO2_JPL-2010(2011)_298K_245-415nm(rec).txt', '1'),
    ('Jn50', 'Cross_Section_ICl_JPL-2010(2011)_298K_210-600nm(rec).txt', '1'),
    ('Jn51', 'Cross_Section_IBr_JPL-2010(2011)_298K_220-600nm(rec).txt', '1'),
    ('Jn52', 'Cross_Section_FURFURAL.csv',      '0.5'),
    ('Jn53', 'Cross_Section_ClNO_IUPAC2006.csv','1'),
    ('Jn54', 'Cross_Section_ClONO_IUPAC2006.csv','1'),
    ('Jn55', 'CH2ClCHO_JPL-2010(2011)_298K_240-357nm(rec).txt', '1'),
    ('Jn56', 'CH3C(O)CH2Cl_JPL-2010(2011)_296K_210-360nm(rec).txt', '1'),
]


# ── 辅助函数 ──

def read_csv(path):
    """读取 CSV/TXT 文件，返回 (wl, [col0, col1, ...])。支持逗号、空格、制表符分隔。"""
    wl, cols = [], []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('%'):
                continue
            # 统一分隔符
            for sep in [',', '\t', ' ']:
                if sep in line:
                    parts = line.split(sep)
                    break
            else:
                continue
            parts = [p.strip() for p in parts if p.strip()]
            if len(parts) < 2:
                continue
            try:
                w = float(parts[0])
                if w <= 0:
                    continue
            except ValueError:
                continue
            vals = []
            for p in parts[1:]:
                try:
                    vals.append(float(p))
                except ValueError:
                    vals.append(float('nan'))
            if not vals:
                continue
            wl.append(w)
            if not cols:
                cols = [[] for _ in vals]
            for j, v in enumerate(vals):
                cols[j].append(v)
    return wl, cols


def write_csv2(path, wl, vals):
    """写入 2 列 CSV。"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        for w, v in zip(wl, vals):
            f.write(f"{w:.8g},{v:.8e}\n")


def handle_generic_2col(path):
    """简单 2 列 CSV 直读。"""
    wl, cols = read_csv(path)
    return wl, cols[0] if cols else []


def handle_generic_3col(path, T1=220, T2=298):
    """3 列 CSV，在 T1、T2 间线性内插到 TARGET_T。"""
    wl, cols = read_csv(path)
    if len(cols) < 2:
        return wl, cols[0] if cols else []
    frac = (TARGET_T - T1) / (T2 - T1) if T2 != T1 else 0
    frac = max(0, min(1, frac))
    vals = []
    for i in range(len(wl)):
        v1 = cols[0][i] if i < len(cols[0]) else float('nan')
        v2 = cols[1][i] if i < len(cols[1]) else float('nan')
        if not math.isnan(v1) and not math.isnan(v2):
            vals.append(v1 + frac * (v2 - v1))
        elif not math.isnan(v2):
            vals.append(v2)
        elif not math.isnan(v1):
            vals.append(v1)
        else:
            vals.append(0)
    return wl, vals


# ── MATLAB 函数句柄的手动实现 (T=298K 预计算) ──

def handle_O3_JPL_CS():
    wl, cols = read_csv(CS_SRC / "Cross_Section_O3_JPL.csv")
    # T1=218K(col1), T2=295K(col2), 目标 298K → 用 col2
    vals = []
    for i in range(len(wl)):
        v2 = cols[1][i] if len(cols) > 1 and i < len(cols[1]) and not math.isnan(cols[1][i]) else float('nan')
        v1 = cols[0][i] if len(cols) > 0 and i < len(cols[0]) and not math.isnan(cols[0][i]) else float('nan')
        if not math.isnan(v2):
            vals.append(v2)
        elif not math.isnan(v1):
            vals.append(v1)
        else:
            vals.append(0)
    return wl, vals


def handle_O3_O1D_QY():
    """Quantum_Yield_O3_O1D_JPL: 分段多项式。简化: 读取 IUPAC CSV。"""
    try:
        return handle_generic_2col(QY_SRC / "Quantum_Yield_O3_O1D_IUPAC.csv")
    except:
        return handle_generic_2col(QY_SRC / "Quantum_Yield_O3_O1D_JPL.m")


def handle_O3_O3P_QY():
    """QY_O3P = 1 - QY_O1D。"""
    wl, qy_o1d = handle_O3_O1D_QY()
    return wl, [1.0 - q for q in qy_o1d]


def handle_NO2_CS():
    wl, cols = read_csv(CS_SRC / "Cross_Section_NO2.csv")
    # 值单位是 10^-20 cm², 需要除以 1e20
    vals = [cols[0][i] / 1e20 if i < len(cols[0]) else 0 for i in range(len(wl))]
    return wl, vals


def handle_NO2_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_NO2.csv")


def handle_H2O2_CS():
    wl, cols = read_csv(CS_SRC / "Cross_Section_H2O2.csv")
    vals = list(cols[0])
    # λ > 260nm 有温度依赖，T=298K 时简化处理
    return wl, vals


def handle_NO3_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_NO3.csv")


def handle_HNO3_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_HNO3.csv")


def handle_HCHO_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_HCHO.csv")


def handle_CH3COCH3_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_CH3COCH3.csv")


def handle_CH3NO3_CS():
    return handle_generic_3col(CS_SRC / "Cross_Section_CH3NO3.csv", 220, 298)


def handle_C2H5NO3_CS():
    return handle_generic_3col(CS_SRC / "Cross_Section_C2H5NO3.csv", 220, 298)


def handle_ClONO2_CS():
    return handle_generic_3col(CS_SRC / "Cross_Section_ClONO2.csv", 220, 296)


def handle_PAN_CS():
    return handle_generic_3col(CS_SRC / "Cross_Section_PAN.csv", 220, 298)


def handle_IC3H7NO3_CS():
    return handle_generic_3col(CS_SRC / "Cross_Section_IC3H7NO3.csv", 220, 298)


def handle_ClNO2_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_ClNO2.csv")


def handle_Cl2_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_Cl2.csv")


def handle_BrONO2_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_BrONO2.csv")


def handle_CHBr3_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_CHBr3.csv")


def handle_Br2_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_Br2.csv")


def handle_N2O5_CS():
    return handle_generic_2col(CS_SRC / "Cross_Section_N2O5.csv")


def handle_GLYOX_H2_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_GLYOX_JPL.csv")


def handle_HCHO_HCO_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_HCHO_H2.m")


def handle_HCHO_H2_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_HCHO_H2.csv")


def handle_CH3CHO_CH3_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_CH3CHO.csv")


def handle_C3H7CHO_C3H7_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_C3H7CHO_C2H4.m")


def handle_C3H7CHO_C2H4_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_C3H7CHO.csv")


def handle_MVK_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_MVK.m")


def handle_GLYOX_HCHO_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_GLYOX_JPL.csv")


def handle_GLYOX_HCO_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_GLYOX_JPL.csv")


def handle_MGLYOX_QY():
    """简化分段函数。"""
    wl, vals = [], []
    for w in range(200, 481):
        wl.append(float(w))
        if w < 290:
            vals.append(0.5)
        elif w < 350:
            vals.append(0.5 * (350 - w) / 60.0)
        elif w < 400:
            vals.append(0.05)
        else:
            vals.append(0.01)
    return wl, vals


def handle_Acrolein_QY():
    return handle_generic_2col(QY_SRC / "Quantum_Yield_CH3CHO.csv")


# ── 处理函数映射表 ──
M_FUNCTION_HANDLERS = {
    'Cross_Section_O3_JPL':          (handle_O3_JPL_CS,      'Cross_Section_O3_JPL_precomp.csv', False),
    'Quantum_Yield_O3_O1D_JPL':      (handle_O3_O1D_QY,      'Quantum_Yield_O3_O1D_JPL_precomp.csv', True),
    'Quantum_Yield_O3_O3P_JPL':      (handle_O3_O3P_QY,      'Quantum_Yield_O3_O3P_JPL_precomp.csv', True),
    'Cross_Section_NO2':             (handle_NO2_CS,         'Cross_Section_NO2_precomp.csv', False),
    'Quantum_Yield_NO2':             (handle_NO2_QY,         'Quantum_Yield_NO2_precomp.csv', True),
    'Cross_Section_H2O2':            (handle_H2O2_CS,        'Cross_Section_H2O2_precomp.csv', False),
    'Cross_Section_NO3':             (handle_NO3_CS,         'Cross_Section_NO3_precomp.csv', False),
    'Cross_Section_HNO3':            (handle_HNO3_CS,        'Cross_Section_HNO3_precomp.csv', False),
    'Cross_Section_HCHO':            (handle_HCHO_CS,        'Cross_Section_HCHO_precomp.csv', False),
    'Cross_Section_CH3COCH3':        (handle_CH3COCH3_CS,    'Cross_Section_CH3COCH3_precomp.csv', False),
    'Cross_Section_CH3NO3':          (handle_CH3NO3_CS,      'Cross_Section_CH3NO3_precomp.csv', False),
    'Cross_Section_C2H5NO3':         (handle_C2H5NO3_CS,     'Cross_Section_C2H5NO3_precomp.csv', False),
    'Cross_Section_ClONO2':          (handle_ClONO2_CS,      'Cross_Section_ClONO2_precomp.csv', False),
    'Cross_Section_PAN':             (handle_PAN_CS,         'Cross_Section_PAN_precomp.csv', False),
    'Cross_Section_IC3H7NO3':        (handle_IC3H7NO3_CS,    'Cross_Section_IC3H7NO3_precomp.csv', False),
    'Cross_Section_ClNO2':           (handle_ClNO2_CS,       'Cross_Section_ClNO2_precomp.csv', False),
    'Cross_Section_Cl2':             (handle_Cl2_CS,         'Cross_Section_Cl2_precomp.csv', False),
    'Cross_Section_BrONO2':          (handle_BrONO2_CS,      'Cross_Section_BrONO2_precomp.csv', False),
    'Cross_Section_CHBr3':           (handle_CHBr3_CS,       'Cross_Section_CHBr3_precomp.csv', False),
    'Cross_Section_Br2':             (handle_Br2_CS,         'Cross_Section_Br2_precomp.csv', False),
    'Cross_Section_N2O5':            (handle_N2O5_CS,        'Cross_Section_N2O5_precomp.csv', False),
    'Quantum_Yield_HCHO_HCO':        (handle_HCHO_HCO_QY,    'Quantum_Yield_HCHO_HCO_precomp.csv', True),
    'Quantum_Yield_HCHO_H2':         (handle_HCHO_H2_QY,     'Quantum_Yield_HCHO_H2_precomp.csv', True),
    'Quantum_Yield_CH3CHO_CH3':      (handle_CH3CHO_CH3_QY,  'Quantum_Yield_CH3CHO_CH3_precomp.csv', True),
    'Quantum_Yield_C3H7CHO_C3H7':    (handle_C3H7CHO_C3H7_QY,'Quantum_Yield_C3H7CHO_C3H7_precomp.csv', True),
    'Quantum_Yield_C3H7CHO_C2H4':    (handle_C3H7CHO_C2H4_QY,'Quantum_Yield_C3H7CHO_C2H4_precomp.csv', True),
    'Quantum_Yield_MVK':             (handle_MVK_QY,         'Quantum_Yield_MVK_precomp.csv', True),
    'Quantum_Yield_GLYOX_H2':        (handle_GLYOX_H2_QY,    'Quantum_Yield_GLYOX_H2_precomp.csv', True),
    'Quantum_Yield_GLYOX_HCHO':      (handle_GLYOX_HCHO_QY,  'Quantum_Yield_GLYOX_HCHO_precomp.csv', True),
    'Quantum_Yield_GLYOX_HCO':       (handle_GLYOX_HCO_QY,   'Quantum_Yield_GLYOX_HCO_precomp.csv', True),
    'Quantum_Yield_MGLYOX':          (handle_MGLYOX_QY,      'Quantum_Yield_MGLYOX_precomp.csv', True),
    'Quantum_Yield_Acrolein':        (handle_Acrolein_QY,    'Quantum_Yield_Acrolein_precomp.csv', True),
    'Quantum_Yield_CH3COCH3_CH3CO':  (None, 'Quantum_Yield_CH3COCH3_CH3CO_scalar.csv', True),
    'Quantum_Yield_CH3COCH3_CO':     (None, 'Quantum_Yield_CH3COCH3_CO_scalar.csv', True),
    'Quantum_Yield_CH3CHO_CH4':      (None, 'Quantum_Yield_CH3CHO_CH4_scalar.csv', True),
    'Quantum_Yield_CH3CHO_CH3CO':    (None, 'Quantum_Yield_CH3CHO_CH3CO_scalar.csv', True),
}


# ── 全局变量（由 main 设置）──
CS_SRC, QY_SRC, CS_DST, QY_DST, MAP_FILE, TARGET_T, TARGET_P = [None]*7


def resolve_source(source_str, is_qy=False):
    """解析 CS 或 QY 来源字符串，返回 (type, filename)。
    type: 0=标量, 1=2列CSV, 3=TXT"""
    s = source_str.strip()

    if s.startswith('@'):
        fname = s[1:]
        if fname in M_FUNCTION_HANDLERS:
            handler, outfile, _is_qy = M_FUNCTION_HANDLERS[fname]
            if handler:
                dst_dir = QY_DST if _is_qy else CS_DST
                dst_path = dst_dir / outfile
                if not dst_path.exists():
                    print(f"  预计算 {fname} → {outfile}...")
                    wl, vals = handler()
                    write_csv2(dst_path, wl, vals)
                return (1, outfile)
            else:
                return (0, outfile)
        else:
            print(f"  警告: 未知函数句柄 @{fname}，跳过")
            return (1, f"UNKNOWN_{fname}")

    elif s.endswith('.csv'):
        # 直接使用
        return (1, s)
    elif s.endswith('.txt'):
        # TXT 文件
        return (3, s)
    else:
        # 标量 QY
        return (0, s)


def main():
    global CS_SRC, QY_SRC, CS_DST, QY_DST, MAP_FILE, TARGET_T, TARGET_P

    import argparse as ap
    parser = ap.ArgumentParser(
        description='预计算 BottomUp 光解的 CS/QY CSV 并生成反应映射文件',
        formatter_class=ap.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s /path/to/F0AM/Chem/Photolysis /path/to/output 298.0
  %(prog)s Chem/Photolysis database/photolysis/bottomup
        """
    )
    parser.add_argument('photolysis_dir', help='F0AM Chem/Photolysis 目录路径')
    parser.add_argument('output_dir', help='输出根目录')
    parser.add_argument('temperature', nargs='?', type=float, default=298.0,
                        help='预计算温度 (K), 默认 298.0')
    args = parser.parse_args()

    TARGET_T = args.temperature
    TARGET_P = 1013.0

    photolysis = Path(args.photolysis_dir)
    if not photolysis.is_dir():
        print(f"错误: 目录不存在 '{photolysis}'", file=sys.stderr)
        sys.exit(1)

    CS_SRC = photolysis / "CrossSections"
    QY_SRC = photolysis / "QuantumYields"
    if not CS_SRC.is_dir():
        print(f"错误: CrossSections 目录不存在 '{CS_SRC}'", file=sys.stderr)
        sys.exit(1)
    if not QY_SRC.is_dir():
        print(f"警告: QuantumYields 目录不存在 '{QY_SRC}'", file=sys.stderr)

    out_dir = Path(args.output_dir)
    CS_DST = out_dir / "CrossSections"
    QY_DST = out_dir / "QuantumYields"
    MAP_FILE = out_dir / "bottomup_jmap.dat"

    CS_DST.mkdir(parents=True, exist_ok=True)
    QY_DST.mkdir(parents=True, exist_ok=True)

    print(f"=== BottomUp J 值预计算 ===")
    print(f"目标: T={TARGET_T}K, P={TARGET_P}mbar")
    print(f"输出目录: CS={CS_DST}, QY={QY_DST}")

    # 处理每个反应
    map_lines = []
    for jname, cs_src, qy_src in J_BOTTOMUP_REACTIONS:
        print(f"处理 {jname}...")
        cs_type, cs_file = resolve_source(cs_src, is_qy=False)
        qy_type, qy_file = resolve_source(qy_src, is_qy=True)
        map_lines.append(f"{jname}\t{cs_file}\t{cs_type}\t{qy_file}\t{qy_type}")

    # 写映射文件
    with open(MAP_FILE, 'w') as f:
        f.write("# BottomUp photolysis reaction mapping\n")
        f.write("# Format: JNAME  CS_FILE  CS_TYPE  QY_FILE  QY_TYPE\n")
        f.write("# CS_TYPE: 1=2-column CSV, 2=3-column CSV(T-interp), 3=TXT\n")
        f.write("# QY_TYPE: 0=scalar, 1=2-column CSV, 2=3-column CSV, 3=TXT\n")
        f.write(f"# Pre-computed at T={TARGET_T}K, P={TARGET_P}mbar\n\n")
        for line in map_lines:
            f.write(line + '\n')

    print(f"\n映射文件已写入: {MAP_FILE}")
    print(f"总反应数: {len(map_lines)}")
    print("完成!")


if __name__ == '__main__':
    main()
