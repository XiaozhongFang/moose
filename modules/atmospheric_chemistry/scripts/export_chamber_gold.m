% export_chamber_gold.m
% Run all chamber scenarios and export/copy gold CSV files for MOOSE CSVDiff tests.
%
% Usage:
%   cd('D:/BaiduSyncdisk/Code/F0AM/Scripts/Tutorials/ExampleSetup_Chamber')
%   addpath(genpath('D:/BaiduSyncdisk/Code/F0AM'))
%   export_chamber_gold

function export_chamber_gold()
    f0am_root = 'D:/BaiduSyncdisk/Code/F0AM';
    runs_dir = fullfile(f0am_root, 'Runs');
    moose_gold_dir = 'W:/home/fangxiaozhong/git_repo/moose/modules/atmospheric_chemistry/test/tests/actions/gold';

    addpath(genpath(f0am_root));

    target_times = 1000:1000:10000;

    % Main run: S1/S2/S3
    Met = {'P',1013; 'T',298; 'RH',10; 'LFlux','ExampleLightFlux.txt'; 'jcorr',1};
    InitConc = {'C5H8',10,0; 'NO2',[0.1;1;10],0; 'H2O2',200,0};
    ChemFiles = {'MCMv331_K(Met)'; 'MCMv331_J(Met,1)'; 'MCMv331_Inorg_Isoprene'; 'CH3ONO_hv'};
    BkgdConc = {'DEFAULT',0};

    ModelOptions.Verbose = 1;
    ModelOptions.EndPointsOnly = 0;
    ModelOptions.LinkSteps = 0;
    ModelOptions.IntTime = 10800;
    ModelOptions.SavePath = fullfile(runs_dir, 'AllScenariosOutput.mat');
    ModelOptions.GoParallel = 0;

    fprintf('=== Running F0AM main simulation (S1/S2/S3) ===\n');
    S = F0AM_ModelCore(Met, InitConc, ChemFiles, BkgdConc, ModelOptions);

    air_den = resolve_air_density(S);
    fprintf('  Air density used for ppb->molec/cm^3: %.6e\n', air_den);

    SplitRun(S, 'step');
    if ~exist('S1', 'var') || ~exist('S2', 'var') || ~exist('S3', 'var')
        error('SplitRun did not create S1/S2/S3 in caller workspace.');
    end

    write_scenario_csv(S1, fullfile(runs_dir, 'F0AM_S1_gold.csv'), target_times, air_den, 'S1');
    write_scenario_csv(S2, fullfile(runs_dir, 'F0AM_S2_gold.csv'), target_times, air_den, 'S2');
    write_scenario_csv(S3, fullfile(runs_dir, 'F0AM_S3_gold.csv'), target_times, air_den, 'S3');

    % Restart run: S2b (jcorr=10)
    fprintf('=== Running F0AM restart simulation (S2b) ===\n');
    [InitConc2, Met2] = Run2Init(S2, length(S2.Time));
    j_idx = ismember(Met2(:,1), 'jcorr');
    Met2{j_idx,2} = 10;

    ModelOptions2 = ModelOptions;
    ModelOptions2.IntTime = 3600;
    ModelOptions2.SavePath = fullfile(runs_dir, 'AllScenarios_S2b.mat');

    S2b = F0AM_ModelCore(Met2, InitConc2, ChemFiles, BkgdConc, ModelOptions2);

    t_offset = S2.Time(end);
    S2b.Time = S2b.Time + t_offset;

    s2b_duration = S2b.Time(end) - t_offset;
    target_s2b = 1000:1000:min(10000, floor(s2b_duration));
    target_s2b_abs = t_offset + target_s2b;

    write_scenario_csv(S2b, fullfile(runs_dir, 'F0AM_S2b_gold.csv'), target_s2b_abs, air_den, 'S2b');

    if ~isfolder(moose_gold_dir)
        mkdir(moose_gold_dir);
    end

    src_names = {
        'F0AM_S1_gold.csv'
        'F0AM_S2_gold.csv'
        'F0AM_S3_gold.csv'
        'F0AM_S2b_gold.csv'
    };

    dst_names = {
        'vs_F0AM_chamber_S1_box.csv'
        'vs_F0AM_chamber_S2_box.csv'
        'vs_F0AM_chamber_S3_box.csv'
        'vs_F0AM_chamber_S2b_box.csv'
    };

    for i = 1:numel(src_names)
        src = fullfile(runs_dir, src_names{i});
        dst = fullfile(moose_gold_dir, dst_names{i});

        if exist(src, 'file') ~= 2
            error('Missing generated file: %s', src);
        end

        copyfile(src, dst);
        fprintf('Copied %s -> %s\n', src_names{i}, dst);
    end

    fprintf('\nChamber gold CSV update complete.\n');
end

function air_den = resolve_air_density(S)
    air_den = NaN;

    % Case 1: struct style, e.g. S.Met.M
    if isfield(S, 'Met') && isstruct(S.Met) && isfield(S.Met, 'M')
        m = S.Met.M;
        if isnumeric(m) && ~isempty(m)
            air_den = m(1);
            return;
        end
    end

    % Case 2: cell style, e.g. S.Met = {'M', value; ...}
    if isfield(S, 'Met') && iscell(S.Met) && size(S.Met,2) >= 2
        idx = find(strcmp(S.Met(:,1), 'M'), 1);
        if ~isempty(idx)
            v = S.Met{idx,2};
            if isnumeric(v) && ~isempty(v)
                air_den = v(1);
                return;
            end
        end
    end

    % Safe fallback for this chamber setup (T=298K, P=1013mbar)
    air_den = 2.46e19;
    warning('Met.M not found. Fallback air density %.3e is used.', air_den);
end

function write_scenario_csv(Sstruct, filename, target_times, air_den, tag)
    Cnames = Sstruct.Cnames;
    t = Sstruct.Time;

    fid = fopen(filename, 'w');
    if fid < 0
        error('Cannot open file for write: %s', filename);
    end

    fprintf(fid, 'time');
    for si = 1:length(Cnames)
        fprintf(fid, ',%s', Cnames{si});
    end
    fprintf(fid, '\n');

    for ti = 1:length(target_times)
        tt = target_times(ti);
        [~, idx] = min(abs(t - tt));
        actual_t = t(idx);
        fprintf(fid, '%.0f', actual_t);
        for si = 1:length(Cnames)
            v = Sstruct.Conc.(Cnames{si});
            v_molec = v(idx) * air_den / 1e9;
            fprintf(fid, ',%.12e', v_molec);
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
    fprintf('Wrote %s CSV: %s\n', tag, filename);
end
