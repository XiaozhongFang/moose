% export_f0am_all_scenarios.m
% 完整复现 ExampleSetup_Chamber.m 的所有 4 个场景并导出 gold CSV
%
% 场景:
%   S1:  NO2=0.1ppb, jcorr=1   (low-NOx)
%   S2:  NO2=1ppb,   jcorr=1   (mid-NOx)
%   S3:  NO2=10ppb,  jcorr=1   (hi-NOx)
%   S2b: NO2=1ppb,   jcorr=10  (mid-NOx, 10×光强, 从S2末态重启)
%
% 用法:
%   cd('D:\BaiduSyncdisk\Code\F0AM');
%   addpath(genpath(pwd));
%   export_f0am_all_scenarios
%
% 输出:
%   Runs/F0AM_S{1,2,3}_gold.csv  — 主运行 3 个步骤 (molec/cm³)
%   Runs/F0AM_S2b_gold.csv       — EVENTS 续算 (molec/cm³)

function export_f0am_all_scenarios()
    addpath(genpath('D:\BaiduSyncdisk\Code\F0AM'));
    
    % S1/S2/S3 的目标输出时间点 (相对于各步起始时间为 0~3h)
    target_times = 1000:1000:10000;
    outdir = 'D:\BaiduSyncdisk\Code\F0AM\Runs\';
    
    % === Main run: 3 个 NO2 水平 ===
    Met = {'P',1013; 'T',298; 'RH',10; 'LFlux','ExampleLightFlux.txt'; 'jcorr',1};
    InitConc = {'C5H8',10,0; 'NO2',[0.1;1;10],0; 'H2O2',200,0};
    ChemFiles = {'MCMv331_K(Met)'; 'MCMv331_J(Met,1)'; 'MCMv331_Inorg_Isoprene'; 'CH3ONO_hv'};
    BkgdConc = {'DEFAULT',0};
    ModelOptions.Verbose = 1;
    ModelOptions.EndPointsOnly = 0;
    ModelOptions.LinkSteps = 0;
    ModelOptions.IntTime = 10800;
    ModelOptions.SavePath = [outdir 'AllScenariosOutput.mat'];
    ModelOptions.GoParallel = 0;
    
    fprintf('=== Running F0AM main simulation (S1/S2/S3) ===\n');
    S = F0AM_ModelCore(Met, InitConc, ChemFiles, BkgdConc, ModelOptions);
    
    % 从 Met 读取实际空气数密度 (随 P,T 动态计算)
    M_idx = find(strcmp(S.Met(:,1), 'M'));
    if isempty(M_idx)
        error('Met.M not found in output. Cannot determine air density.');
    end
    air_den = S.Met{M_idx, 2};  % molecules/cm³
    fprintf('  Air density: %.4e molec/cm³\n', air_den);
    
    % SplitRun 拆分为 S1/S2/S3
    S_split = SplitRun(S, 'step');
    S1 = S_split{1};
    S2 = S_split{2};
    S3 = S_split{3};
    
    % 导出 S1/S2/S3
    export_scenario(S1, [outdir 'F0AM_S1_gold.csv'], target_times, air_den, 'S1 (NO2=0.1ppb)');
    export_scenario(S2, [outdir 'F0AM_S2_gold.csv'], target_times, air_den, 'S2 (NO2=1ppb)');
    export_scenario(S3, [outdir 'F0AM_S3_gold.csv'], target_times, air_den, 'S3 (NO2=10ppb)');
    
    % === EVENTS: S2b (重启, jcorr=10) ===
    fprintf('\n=== Running EVENTS (S2b): jcorr=10, restart from S2 end ===\n');
    [InitConc2, Met2] = Run2Init(S2, length(S2.Time));
    loc = ismember(Met2(:,1), 'jcorr');
    Met2{loc,2} = 10;
    
    ModelOptions2 = ModelOptions;
    ModelOptions2.IntTime = 3600;  % 续跑 1h
    ModelOptions2.SavePath = [outdir 'AllScenarios_S2b.mat'];
    
    S2b = F0AM_ModelCore(Met2, InitConc2, ChemFiles, BkgdConc, ModelOptions2);
    
    % S2b 时间平移: 内部从 t=0 开始 → 对齐到 S2 末态
    t_offset = S2.Time(end);  % ≈10800s
    S2b.Time = S2b.Time + t_offset;
    
    % S2b 的目标时间点: 相对于重启时刻的偏移量
    % 取 0~S2b.IntTime 范围内的标准间隔
    s2b_duration = S2b.Time(end) - t_offset;  % ≈3600s
    target_s2b = 1000:1000:min(10000, s2b_duration);
    % 映射到绝对时间
    target_s2b_abs = t_offset + target_s2b;
    
    export_scenario(S2b, [outdir 'F0AM_S2b_gold.csv'], target_s2b_abs, air_den, ...
        sprintf('S2b (NO2=1ppb, jcorr=10, t_offset=%.0fs)', t_offset));
    
    save([outdir 'AllScenariosOutput.mat'], 'S1', 'S2', 'S3', 'S2b', 'S');
    fprintf('\n=== ALL DONE ===\n');
end

function export_scenario(Sstruct, filename, target_times, air_den, label)
    Cnames = Sstruct.Cnames;
    t = Sstruct.Time;
    
    % 验证目标时间是否在数据范围内
    if target_times(end) > t(end)
        warning('export_scenario: %s — target time %.0f > last data time %.0f', ...
            label, target_times(end), t(end));
    end
    if target_times(1) < t(1)
        warning('export_scenario: %s — target time %.0f < first data time %.0f', ...
            label, target_times(1), t(1));
    end
    
    fid = fopen(filename, 'w');
    fprintf(fid, 'time');
    for si = 1:length(Cnames)
        fprintf(fid, ',%s', Cnames{si});
    end
    fprintf(fid, '\n');
    
    for ti = 1:length(target_times)
        tt = target_times(ti);
        [~, idx] = min(abs(t - tt));
        % 输出实际匹配的时间（用于验证）
        actual_t = t(idx);
        fprintf(fid, '%.0f', actual_t);
        for si = 1:length(Cnames)
            val = Sstruct.Conc.(Cnames{si});
            val_molec = val(idx) * air_den / 1e9;  % ppb → molec/cm³
            fprintf(fid, ',%.12e', val_molec);
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('  %s → %s (%d sp × %d tp)\n', label, filename, length(Cnames), length(target_times));
end
