function run_hmm_on_population(config_path, popcurrent_path, output_path) % Optional
    : open parallel pool if not already open if isempty (gcp('nocreate')) parpool('Processes', 4);
end

    % Setup paths homedir = getuserdir();
addpath(genpath(fullfile('..', 'HMM-MAR')));
addpath(genpath(fullfile('..', 'HMM-MAR', 'HMM-MAR')));
addpath(genpath(fullfile('..', 'nutmeg')));

% Load config str = fileread(config_path);
data = jsondecode(str);

% Load data load(popcurrent_path, 'popcurrent');

sizePop = data.GA_config.size_pop;
sampling_freq = data.solve_options.sampling_rate;
numTrials = data.solve_options.mc_trials;

hmmstats = NaN(4, 2, sizePop);

% Optionally save intermediate stats save(fullfile(output_path, 'stats.mat'), 'hmmstats');

disp('Running HMM on all phenotypes in this generation');

parfor i = 1 : sizePop fprintf('Running HMM on phenotype %d\n', i);
hmmstats( :, :, i) = supress_gather_output(i, popcurrent, numTrials, sampling_freq);
end

    stats = hmmstats;
save(fullfile(output_path, 'HMMStats.mat'), 'stats');
end

    function stats = supress_gather_output(i, popcurrent, numTrials, sampling_freq) try[~, stats] =
        evalc('hmm_burst_detect_and_stats(squeeze(popcurrent(:,:,i)), numTrials, sampling_freq);');
catch stats = NaN(4, 2);
end end