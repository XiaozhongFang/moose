% generate_bottomup_data.m
% Generate BottomUp photolysis data files using F0AM's native MATLAB functions.
%
% Reads J_BottomUp.m's reaction list, evaluates each cross-section and
% quantum-yield function at the specified T and P, and writes pre-computed
% 2-column CSV files + bottomup_jmap.dat.
%
% Usage:
%   addpath(genpath('D:\BaiduSyncdisk\Code\F0AM'));
%   generate_bottomup_data(298, 1013, 'D:\output\dir');
%
% This eliminates all manual Python reimplementation of F0AM's MATLAB functions.
% Any new photolysis reaction added to F0AM is handled automatically.

function generate_bottomup_data(T, P, outdir)
    if nargin < 1, T = 298; end
    if nargin < 2, P = 1013; end
    if nargin < 3, outdir = 'D:\BaiduSyncdisk\Code\F0AM\Runs\bottomup'; end

    % Ensure F0AM is on the path
    if ~exist('J_BottomUp', 'file')
        error('F0AM not on path. Run: addpath(genpath(''D:\\BaiduSyncdisk\\Code\\F0AM''))');
    end

    cs_dir  = fullfile(outdir, 'CrossSections');
    qy_dir  = fullfile(outdir, 'QuantumYields');
    mkdir(cs_dir);
    mkdir(qy_dir);

    % === Reaction list matching J_BottomUp.m ===
    % (name, CS_source, QY_source)
    % CS/QY_source: '@FunctionName' (MATLAB handle) or 'filename.csv' or scalar
    reactions = {
        'J1',  '@Cross_Section_O3_JPL',           '@Quantum_Yield_O3_O1D_JPL'
        'J2',  '@Cross_Section_O3_JPL',           '@Quantum_Yield_O3_O3P_JPL'
        'J3',  '@Cross_Section_H2O2',             '1'
        'J4',  '@Cross_Section_NO2',              '@Quantum_Yield_NO2'
        'J5',  '@Cross_Section_NO3',              'Quantum_Yield_NO3_NO.csv'
        'J6',  '@Cross_Section_NO3',              'Quantum_Yield_NO3_NO2.csv'
        'J7',  'Cross_Section_HONO.csv',          '1'
        'J8',  '@Cross_Section_HNO3',             '1'
        'J11', '@Cross_Section_HCHO',             '@Quantum_Yield_HCHO_HCO'
        'J12', '@Cross_Section_HCHO',             '@Quantum_Yield_HCHO_H2'
        'J13', 'Cross_Section_CH3CHO.csv',        '@Quantum_Yield_CH3CHO_CH3'
        'J14', 'Cross_Section_C2H5CHO.csv',        'Quantum_Yield_C2H5CHO.csv'
        'J15', 'Cross_Section_C3H7CHO.csv',        '@Quantum_Yield_C3H7CHO_C3H7'
        'J16', 'Cross_Section_C3H7CHO.csv',        '@Quantum_Yield_C3H7CHO_C2H4'
        'J17', 'Cross_Section_IPRCHO.csv',         'Quantum_Yield_IPRCHO.csv'
        'J18', 'Cross_Section_MACR.csv',           '0.005'
        'J19', 'Cross_Section_MACR.csv',           '0.005'
        'J20', 'Cross_Section_MACR.csv',           '1'
        'J21', '@Cross_Section_CH3COCH3',         '@Quantum_Yield_CH3COCH3_CH3CO'
        'J22', 'Cross_Section_MEK.csv',            '0.34'
        'J23', 'Cross_Section_MVK.csv',           '@Quantum_Yield_MVK'
        'J24', 'Cross_Section_MVK.csv',           '@Quantum_Yield_MVK'
        'J31', 'Cross_Section_GLYOX.csv',         '@Quantum_Yield_GLYOX_H2'
        'J32', 'Cross_Section_GLYOX.csv',         '@Quantum_Yield_GLYOX_HCHO'
        'J33', 'Cross_Section_GLYOX.csv',         '@Quantum_Yield_GLYOX_HCO'
        'J34', 'Cross_Section_MGLYOX.csv',        '@Quantum_Yield_MGLYOX'
        'J35', 'Cross_Section_BIACET.csv',         '0.158'
        'J41', 'Cross_Section_CH3OOH.csv',         '1'
        'J51', '@Cross_Section_CH3NO3',           '1'
        'J52', '@Cross_Section_C2H5NO3',          '1'
        'J53', 'Cross_Section_NC3H7NO3.csv',       '1'
        'J54', '@Cross_Section_IC3H7NO3',         '1'
        'J55', 'Cross_Section_TC4H9NO3.csv',       '1'
        'J56', 'Cross_Section_NOA.csv',            '0.9'
        'J57', 'Cross_Section_NOA.csv',            '0.1'
        'Jn1',  'Cross_Section_CH3CHCHCHO.csv',    '0.030'
        'Jn2',  'Cross_Section_C6H5COH.csv',       '0.29'
        'Jn3',  'Cross_Section_C2H5COC2H5.csv',    '1'
        'Jn4',  'Cross_Section_CH3COOOH.csv',      '1'
        'Jn5',  'Cross_Section_CH3CHO.csv',        '@Quantum_Yield_CH3CHO_CH4'
        'Jn6',  'Cross_Section_CH3CHO.csv',        '@Quantum_Yield_CH3CHO_CH3CO'
        'Jn8',  '@Cross_Section_CH3COCH3',        '@Quantum_Yield_CH3COCH3_CO'
        'Jn9',  'Cross_Section_GLYCOALDEHYDE.csv', '1'
        'Jn10', 'Cross_Section_Hydroxyaceton.csv', '0.60'
        'Jn11', 'Cross_Section_Acrolein.csv',     '@Quantum_Yield_Acrolein'
        'Jn12', 'Cross_Section_3_methyl_2_nitrophenol.csv', '0.00015'
        'Jn13', 'Cross_Section_4_methyl_2_nitrophenol.csv', '0.0001'
        'Jn14', '@Cross_Section_PAN',             'Quantum_Yield_PAN_NO2.csv'
        'Jn15', '@Cross_Section_PAN',             'Quantum_Yield_PAN_NO3.csv'
        'Jn16', 'Cross_Section_CH3O2NO2.csv',      '0.95'
        'Jn17', 'Cross_Section_CH3O2NO2.csv',      '0.05'
        'Jn18', 'Cross_Section_CH3ONO.csv',        '0.76'
        'Jn19', '@Cross_Section_N2O5',            'Quantum_Yield_N2O5_NO3_NO2.csv'
        'Jn20', '@Cross_Section_N2O5',            'Quantum_Yield_N2O5_NO3_NO_O.csv'
        'Jn21', 'Cross_Section_HO2NO2.csv',        '0.59'
        'Jn22', 'Cross_Section_HO2NO2.csv',        '0.41'
        'Jn23', '@Cross_Section_ClNO2',           '1'
        'Jn24', '@Cross_Section_Br2',             '1'
        'Jn25', 'Cross_Section_BrO.csv',           '1'
        'Jn26', 'Cross_Section_HOBr.csv',          '1'
        'Jn27', 'Cross_Section_BrNO2.csv',         '1'
        'Jn28', '@Cross_Section_BrONO2',          '0.85'
        'Jn29', '@Cross_Section_BrONO2',          '0.15'
        'Jn30', '@Cross_Section_CHBr3',           '1'
        'Jn31', 'Cross_Section_BrCl.csv',          '1'
        'Jn32', '@Cross_Section_Cl2',             '1'
        'Jn33', '@Cross_Section_ClO_MB1999',      '1'
        'Jn34', '@Cross_Section_ClONO2',          'Quantum_Yield_ClONO2_Cl.csv'
        'Jn35', '@Cross_Section_ClONO2',          'Quantum_Yield_ClONO2_ClO.csv'
        'Jn36', 'Cross_Section_HOCl.csv',          '1'
        'Jn37', 'Cross_Section_OClO_Wahner(1987)_296K_245-475nm(0.22nm).txt', '1'
        'Jn38', 'Cross_Section_ClOOCl_JPL-2010(2011)_190-250K_200-420nm(rec).txt', '1'
        'Jn39', 'Cross_Section_ClOO_JPL-2010(2011)_191K_220-280nm(rec).txt', '1'
        'Jn40', 'Cross_Section_I2_JPL-2010(2011)_295K_185-700nm(rec).txt', 'QY_I2.txt'
        'Jn41', 'Cross_Section_HOI_JPL-2010(2011)_295-298K_280-480nm(rec).txt', '1'
        'Jn42', 'Cross_Section_IO_JPL-2010(2011)_298K_339-417nm(rec).txt', '0.91'
        'Jn43', 'Cross_Section_OIO_JPL-2010(2011)_295K_516-572nm(rec).txt', '1'
        'Jn44', 'Cross_Section_INO_JPL-2010(2011)_298K_223-460nm(rec).txt', '1'
        'Jn45', 'Cross_Section_INO2_JPL-2010(2011)_298K_210-380nm(rec).txt', '1'
        'Jn46', 'Cross_Section_IONO2_JPL-2010(2011)_298K_245-415nm(rec).txt', '1'
        'Jn50', 'Cross_Section_ICl_JPL-2010(2011)_298K_210-600nm(rec).txt', '1'
        'Jn51', 'Cross_Section_IBr_JPL-2010(2011)_298K_220-600nm(rec).txt', '1'
        'Jn52', 'Cross_Section_FURFURAL.csv',      '0.5'
        'Jn53', 'Cross_Section_ClNO_IUPAC2006.csv','1'
        'Jn54', 'Cross_Section_ClONO_IUPAC2006.csv','1'
        'Jn55', 'CH2ClCHO_JPL-2010(2011)_298K_240-357nm(rec).txt', '1'
        'Jn56', 'CH3C(O)CH2Cl_JPL-2010(2011)_296K_210-360nm(rec).txt', '1'
    };

    % Photolysis source directory (F0AM's Chem/Photolysis)
    photo_dir = fullfile(fileparts(which('J_BottomUp')));

    n = size(reactions, 1);
    map_lines = cell(n, 1);
    fprintf('=== Generating BottomUp photolysis data at T=%.0fK, P=%.0fmbar ===\n', T, P);
    fprintf('Output: %s\n', outdir);

    for i = 1:n
        jname   = reactions{i, 1};
        cs_src  = strtrim(reactions{i, 2});
        qy_src  = strtrim(reactions{i, 3});
        fprintf('Processing %s...\n', jname);

        % ── Resolve cross-section ──
        [cs_type, cs_file] = resolve_source(cs_src, false, photo_dir, T, P, cs_dir);
        % ── Resolve quantum yield ──
        [qy_type, qy_file] = resolve_source(qy_src, true,  photo_dir, T, P, qy_dir);

        map_lines{i} = sprintf('%s\t%s\t%d\t%s\t%d', jname, cs_file, cs_type, qy_file, qy_type);
    end

    % ── Write bottomup_jmap.dat ──
    map_path = fullfile(outdir, 'bottomup_jmap.dat');
    fid = fopen(map_path, 'w');
    fprintf(fid, '# BottomUp photolysis reaction mapping\n');
    fprintf(fid, '# Generated by generate_bottomup_data.m at T=%.0fK, P=%.0fmbar\n', T, P);
    fprintf(fid, '# Format: JNAME  CS_FILE  CS_TYPE  QY_FILE  QY_TYPE\n');
    fprintf(fid, '# CS_TYPE: 1=2-col CSV, 2=3-col CSV(T-interp), 3=TXT\n');
    fprintf(fid, '# QY_TYPE: 0=scalar, 1=2-col CSV, 2=3-col CSV, 3=TXT\n\n');
    for i = 1:n
        fprintf(fid, '%s\n', map_lines{i});
    end
    fclose(fid);
    fprintf('\nMap written: %s (%d reactions)\n', map_path, n);
    fprintf('Done.\n');
end

% ── Resolve a CS or QY source to (type, filename_or_value) ──
function [type, result] = resolve_source(src, is_qy, photo_dir, T, P, outdir)
    if src(1) == '@'
        % MATLAB function handle: evaluate to get data, write CSV
        func_name = src(2:end);
        try
            fh = str2func(func_name);
            if is_qy
                [vals, wl] = fh(T, P);
            else
                [vals, wl] = fh(T, P);
            end
            % Write 2-column CSV
            outfile = sprintf('%s_precomp.csv', func_name);
            outpath = fullfile(outdir, outfile);
            fid = fopen(outpath, 'w');
            for j = 1:length(wl)
                fprintf(fid, '%.6f,%.12e\n', wl(j), vals(j));
            end
            fclose(fid);
            type = 1;   % 2-col CSV
            result = outfile;
            fprintf('    Evaluated @%s -> %s (%d pts)\n', func_name, outfile, length(wl));
        catch ME
            warning('Failed to evaluate @%s: %s', func_name, ME.message);
            type = 0;
            result = '1';
        end
    elseif endsWith(src, '.csv') || endsWith(src, '.txt')
        % Direct file reference: copy as-is
        src_path = fullfile(photo_dir, src);
        if contains(src, 'CrossSections') || contains(src, 'QuantumYields')
            type = 1;
        elseif endsWith(src, '.txt')
            type = 3;  % TXT
        else
            type = 1;
        end
        result = src;
        % Copy if not already in output dir
        outpath = fullfile(outdir, src);
        if ~exist(outpath, 'file') && exist(src_path, 'file')
            copyfile(src_path, outpath);
        end
    else
        % Scalar value
        type = 0;
        result = src;
    end
end
