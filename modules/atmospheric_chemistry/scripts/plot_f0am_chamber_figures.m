% plot_f0am_chamber_figures.m
% Reproduce the figures produced by F0AM ExampleSetup_Chamber.m and record
% the F0AM compute time for the same four chamber work items used by the MOOSE
% tests: S1, S2, S2b restart, and S3.
%
% Usage from the MOOSE repository root:
%   matlab -batch "addpath('modules/atmospheric_chemistry/scripts'); plot_f0am_chamber_figures"
%
% Optional:
%   plot_f0am_chamber_figures('/path/to/F0AM', 'out/dir')

function timing = plot_f0am_chamber_figures(f0am_root, outdir)
    if nargin < 1 || isempty(f0am_root)
        f0am_root = resolve_f0am_root();
    end
    if nargin < 2 || isempty(outdir)
        outdir = default_output_dir();
    end

    if exist(f0am_root, 'dir') ~= 7
        error('F0AM root not found: %s', f0am_root);
    end
    if exist(outdir, 'dir') ~= 7
        mkdir(outdir);
    end

    addpath(genpath(f0am_root));
    set(0, 'DefaultFigureVisible', 'off');

    Met = {'P',1013; 'T',298; 'RH',10; 'LFlux','ExampleLightFlux.txt'; 'jcorr',1};
    InitConc = {'C5H8',10,0; 'NO2',[0.1;1;10],0; 'H2O2',200,0};
    ChemFiles = {'MCMv331_K(Met)'; 'MCMv331_J(Met,1)'; 'MCMv331_Inorg_Isoprene'; 'CH3ONO_hv'};
    BkgdConc = {'DEFAULT',0};

    runs_dir = fullfile(f0am_root, 'Runs');
    if exist(runs_dir, 'dir') ~= 7
        mkdir(runs_dir);
    end

    ModelOptions.Verbose = 1;
    ModelOptions.EndPointsOnly = 0;
    ModelOptions.LinkSteps = 0;
    ModelOptions.IntTime = 3 * 3600;
    ModelOptions.SavePath = fullfile(runs_dir, 'AllScenariosOutput.mat');
    ModelOptions.GoParallel = 0;

    fprintf('=== Running F0AM main chamber simulation: S1/S2/S3 ===\n');
    tic;
    S = F0AM_ModelCore(Met, InitConc, ChemFiles, BkgdConc, ModelOptions);
    main_seconds = toc;

    SplitRun(S, 'step');
    if ~exist('S1', 'var') || ~exist('S2', 'var') || ~exist('S3', 'var')
        error('SplitRun did not create S1/S2/S3 in this workspace.');
    end

    lnames = {'S1: NO2 = 0.1 ppb', 'S2: NO2 = 1 ppb', 'S3: NO2 = 10 ppb'};
    Splot = {S1, S2, S3};

    PlotConc('C5H8', Splot, 'lnames', lnames);
    save_current_figure(outdir, '01_conc_C5H8');

    PlotConc('OH', Splot, 'unit', 'percc', 'scale', 1e-6, 'lnames', lnames);
    save_current_figure(outdir, '02_conc_OH');

    PlotConc('NO+NO2', Splot, 'unit', 'pptv', 'lnames', lnames);
    save_current_figure(outdir, '03_conc_NOx');

    PlotConcGroup(S3.Cnames(S3.iRO2), S3, 5, 'ptype', 'fill', 'unit', 'ppb', 'name', 'RO_2');
    save_current_figure(outdir, '04_ro2_group_S3');

    PlotRates('PAN', S1, 4, 'ptype', 'fill', 'unit', 'ppt_h', 'sumEq', 1);
    save_current_figure(outdir, '05_rates_PAN_S1');

    pts2avg = S1.Time > 1800 & S1.Time < 3600;
    PlotRatesAvg('HCHO', S1, 5, 'ptype', 'hbar', 'unit', 'ppb_h', 'pts2avg', pts2avg);
    save_current_figure(outdir, '06_rates_avg_HCHO_S1');

    Reactants = {'C5H8', 'NO2', 'H2O2', 'CO', 'CH4'};
    PlotReactivity('OH', S3, Reactants, 'ptype', 'line');
    save_current_figure(outdir, '07_reactivity_OH_S3');

    yieldWindow = [500 1000];
    PlotYield(S1, 'C5H8', {'C5HPALD1', 'C5HPALD2'}, 'twindow', yieldWindow);
    save_current_figure(outdir, '08_yield_C5HPALD_S1');

    fprintf('=== Running F0AM restart chamber simulation: S2b ===\n');
    [InitConc2, Met2] = Run2Init(S2, length(S2.Time));
    j_idx = ismember(Met2(:,1), 'jcorr');
    Met2{j_idx,2} = 10;

    ModelOptions2 = ModelOptions;
    ModelOptions2.IntTime = 3600;
    ModelOptions2.SavePath = fullfile(runs_dir, 'AllScenarios_S2b.mat');

    tic;
    S2b = F0AM_ModelCore(Met2, InitConc2, ChemFiles, BkgdConc, ModelOptions2);
    s2b_seconds = toc;
    S2b.Time = S2b.Time + S2.Time(end);

    PlotConc('C5H8', {S2, S2b});
    save_current_figure(outdir, '09_conc_C5H8_S2_restart');

    timing.main_S1_S2_S3_seconds = main_seconds;
    timing.restart_S2b_seconds = s2b_seconds;
    timing.total_four_tests_seconds = main_seconds + s2b_seconds;

    write_timing_csv(outdir, timing);
    fprintf('F0AM chamber compute total for S1+S2+S2b+S3: %.6f s\n', ...
            timing.total_four_tests_seconds);
    fprintf('Saved figures and timing CSV under: %s\n', outdir);
end

function root = resolve_f0am_root()
    env_root = getenv('F0AM_ROOT');
    if ~isempty(env_root)
        root = env_root;
        return;
    end

    script_dir = fileparts(mfilename('fullpath'));
    candidate = fullfile(script_dir, '..', '..', '..', '.reasonix', 'docs', ...
                         'BaiduSyncdisk', 'Code', 'F0AM');
    if exist(candidate, 'dir') == 7
        root = candidate;
        return;
    end

    error('Set F0AM_ROOT or pass f0am_root explicitly.');
end

function outdir = default_output_dir()
    script_dir = fileparts(mfilename('fullpath'));
    outdir = fullfile(script_dir, '..', 'test', 'tests', 'chamber', 'f0am_figures');
end

function save_current_figure(outdir, name)
    fig = gcf;
    png_path = fullfile(outdir, [name '.png']);
    fig_path = fullfile(outdir, [name '.fig']);
    saveas(fig, png_path);
    savefig(fig, fig_path);
    fprintf('Saved %s\n', png_path);
end

function write_timing_csv(outdir, timing)
    path = fullfile(outdir, 'f0am_chamber_timing.csv');
    fid = fopen(path, 'w');
    if fid < 0
        error('Cannot open timing CSV for write: %s', path);
    end
    fprintf(fid, 'case,seconds\n');
    fprintf(fid, 'main_S1_S2_S3,%.6f\n', timing.main_S1_S2_S3_seconds);
    fprintf(fid, 'restart_S2b,%.6f\n', timing.restart_S2b_seconds);
    fprintf(fid, 'total_current_four_chamber_tests,%.6f\n', timing.total_four_tests_seconds);
    fclose(fid);
end
