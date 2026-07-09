% plot_chamber_f0am_style.m
% Plot chamber CSV outputs with F0AM's own PlotConc/PlotConcGroup helpers.
%
% Usage from the MOOSE repository root:
%   matlab -batch "addpath('modules/atmospheric_chemistry/scripts'); plot_chamber_f0am_style"
%
% Plot a solver benchmark run directory:
%   matlab -batch "addpath('modules/atmospheric_chemistry/scripts'); plot_chamber_f0am_style([], 'modules/atmospheric_chemistry/test/tests/chamber/kpp_chamber/solver_runs/f0am_accuracy_timing')"
%
% Optional full form:
%   plot_chamber_f0am_style('/path/to/F0AM', '/path/to/csv_dir', '/path/to/outdir', {'kpp_rosenbrock','petsc_bdf'})
%
% This script reproduces the concentration-style F0AM chamber plots from CSV
% outputs. F0AM rate/yield plots require per-reaction diagnostics that are not
% present in ordinary chamber CSV files.

function plot_chamber_f0am_style(f0am_root, csv_dir, outdir, solver_names)
    if nargin < 1 || isempty(f0am_root)
        f0am_root = resolve_f0am_root();
    end
    if nargin < 2 || isempty(csv_dir)
        csv_dir = default_chamber_dir();
    end
    if nargin < 3 || isempty(outdir)
        outdir = fullfile(csv_dir, 'f0am_style_figures');
    end
    if nargin < 4
        solver_names = {};
    end

    if exist(f0am_root, 'dir') ~= 7
        error('F0AM root not found: %s', f0am_root);
    end
    if exist(csv_dir, 'dir') ~= 7
        error('CSV directory not found: %s', csv_dir);
    end
    if exist(outdir, 'dir') ~= 7
        mkdir(outdir);
    end

    addpath(genpath(f0am_root));
    set(0, 'DefaultFigureVisible', 'off');

    air_den = 2.4622e19;
    ro2_species = read_ro2_species(default_chamber_mechanism());
    if isempty(solver_names)
        solver_names = discover_solver_names(csv_dir);
    end

    for si = 1:numel(solver_names)
        solver = solver_names{si};
        fprintf('Plotting F0AM-style chamber figures for %s\n', solver);
        [S1, S2, S3, S2b] = load_solver_set(csv_dir, solver, air_den, ro2_species);
        write_plot_set(S1, S2, S3, S2b, outdir, solver);
    end
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

function chamber_dir = default_chamber_dir()
    script_dir = fileparts(mfilename('fullpath'));
    chamber_dir = fullfile(script_dir, '..', 'test', 'tests', 'chamber');
end

function mech = default_chamber_mechanism()
    script_dir = fileparts(mfilename('fullpath'));
    mech = fullfile(script_dir, '..', 'doc', 'content', 'modules', ...
                    'atmospheric_chemistry', 'database', ...
                    'MCMv331_Inorg_Isoprene.fac');
end

function names = discover_solver_names(csv_dir)
    chamber_file = fullfile(csv_dir, 'vs_F0AM_chamber_S1_box.csv');
    if exist(chamber_file, 'file') == 2
        names = {'chamber'};
        return;
    end

    files = dir(fullfile(csv_dir, '*_S1.csv'));
    names = cell(numel(files), 1);
    for i = 1:numel(files)
        names{i} = regexprep(files(i).name, '_S1\.csv$', '');
    end
    if isempty(names)
        error('No chamber CSV set found in %s', csv_dir);
    end
end

function [S1, S2, S3, S2b] = load_solver_set(csv_dir, solver, air_den, ro2_species)
    S1 = load_scenario_csv(scenario_path(csv_dir, solver, 'S1'), air_den, ro2_species, 1);
    S2 = load_scenario_csv(scenario_path(csv_dir, solver, 'S2'), air_den, ro2_species, 2);
    S3 = load_scenario_csv(scenario_path(csv_dir, solver, 'S3'), air_den, ro2_species, 3);

    s2b_path = scenario_path(csv_dir, solver, 'S2b', false);
    if exist(s2b_path, 'file') == 2
        S2b = load_scenario_csv(s2b_path, air_den, ro2_species, 2);
        S2b = offset_restart_time(S2b, S2.Time(end));
    else
        S2b = [];
    end
end

function S = offset_restart_time(S, restart_time)
    % Chamber restart CSVs are stored with relative time for CSVDiff.
    if ~isempty(S.Time) && max(S.Time) <= restart_time
        S.Time = S.Time + restart_time;
    end
end

function path = scenario_path(csv_dir, solver, scenario, required)
    if nargin < 4
        required = true;
    end
    if strcmp(solver, 'chamber')
        path = fullfile(csv_dir, sprintf('vs_F0AM_chamber_%s_box.csv', scenario));
    else
        path = fullfile(csv_dir, sprintf('%s_%s.csv', solver, scenario));
    end
    if required && exist(path, 'file') ~= 2
        error('Missing chamber CSV: %s', path);
    end
end

function S = load_scenario_csv(path, air_den, ro2_species, step_index)
    opts = detectImportOptions(path, 'FileType', 'text', 'Delimiter', ',');
    opts.VariableNamingRule = 'preserve';
    table_data = readtable(path, opts);
    names = table_data.Properties.VariableNames;

    S = struct();
    S.Time = table_data.(names{1});
    S.StepIndex = step_index * ones(size(S.Time));
    S.Cnames = names(2:end)';
    S.Met = struct('M', air_den);
    S.Conc = struct();
    for i = 2:numel(names)
        species_name = names{i};
        S.Conc.(species_name) = table_data.(species_name) .* 1.0e9 ./ air_den;
    end
    S.iRO2 = find(ismember(S.Cnames, ro2_species));
end

function ro2_species = read_ro2_species(mechanism_path)
    if exist(mechanism_path, 'file') ~= 2
        ro2_species = {};
        return;
    end

    text = fileread(mechanism_path);
    token = regexp(text, 'RO2\s*=\s*([\s\S]*?);', 'tokens', 'once');
    if isempty(token)
        ro2_species = {};
        return;
    end

    raw = regexp(token{1}, '\+', 'split');
    ro2_species = cell(0, 1);
    for i = 1:numel(raw)
        name = strtrim(raw{i});
        name = regexprep(name, '[^A-Za-z0-9_]', '');
        if ~isempty(name)
            ro2_species{end + 1, 1} = name; %#ok<AGROW>
        end
    end
end

function write_plot_set(S1, S2, S3, S2b, outdir, prefix)
    lnames = {'low', 'mid', 'hi'};
    Splot = {S1, S2, S3};

    PlotConc('C5H8', Splot, 'lnames', lnames);
    save_current_figure(outdir, prefix, '01_conc_C5H8');

    PlotConc('OH', Splot, 'unit', 'percc', 'scale', 1e-6, 'lnames', lnames);
    save_current_figure(outdir, prefix, '02_conc_OH');

    PlotConc('NO+NO2', Splot, 'unit', 'pptv', 'lnames', lnames);
    save_current_figure(outdir, prefix, '03_conc_NOx');

    if ~isempty(S3.iRO2)
        PlotConcGroup(S3.Cnames(S3.iRO2), S3, 5, 'ptype', 'fill', ...
                      'unit', 'ppb', 'name', 'RO_2');
        save_current_figure(outdir, prefix, '04_ro2_group_S3');
    end

    if ~isempty(S2b)
        PlotConc('C5H8', {S2, S2b});
        save_current_figure(outdir, prefix, '09_conc_C5H8_S2_restart');
    end
end

function save_current_figure(outdir, prefix, name)
    safe_prefix = regexprep(prefix, '[^A-Za-z0-9_]', '_');
    base = fullfile(outdir, [safe_prefix '_' name]);
    saveas(gcf, [base '.png']);
    savefig(gcf, [base '.fig']);
    close(gcf);
    fprintf('Saved %s.png\n', base);
end
