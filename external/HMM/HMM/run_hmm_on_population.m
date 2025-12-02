function run_hmm_on_population(popcurrent_path, output_path, sampling_freq)
% RUN_HMM_ON_POPULATION  Detect CSVs, run HMM, and save individual + combined results.

fprintf('popcurrent_path: %s\n', popcurrent_path);
fprintf('output_path   : %s\n', output_path);
fprintf('sampling_freq : %g\n', sampling_freq);

% Ensure output directory exists
if ~exist(output_path, 'dir')
    mkdir(output_path);
end

% --- Find candidate CSVs ---
files = dir(fullfile(popcurrent_path, 'meg_data_*.csv'));
if isempty(files)
    warning('No files matching meg_data_*.csv found. Falling back to *.csv');
    files = dir(fullfile(popcurrent_path, '*.csv'));
end
if isempty(files)
    error('No CSV files found in: %s', popcurrent_path);
end

% Sort by name for consistent order
[~, sortIdx] = sort(lower({files.name}));
files = files(sortIdx);

% Store concatenated results
all_stats = [];

for k = 1:numel(files)
    in_path = fullfile(files(k).folder, files(k).name);

    % Try to extract index from "meg_data_<n>.csv"
    tok = regexp(files(k).name, '^meg_data_(\d+)\.csv$', 'tokens', 'once');
    if ~isempty(tok)
        out_name = sprintf('HMMStats_%s.csv', tok{1});
    else
        base = regexprep(files(k).name, '\.csv$', '', 'ignorecase');
        out_name = sprintf('HMMStats_%s.csv', base);
    end
    out_path = fullfile(output_path, out_name);

    try
        % Load CSV
        M = readmatrix(in_path);
        [rows, cols] = size(M);
        num_trials = rows;
        fprintf('FILE %d/%d: %s -> %s\n', k, numel(files), files(k).name, out_name);
        fprintf('  rows × cols: %d × %d | num_trials: %d\n', rows, cols, num_trials);

        % Run analysis
        hmm_stats = hmm_burst_detect_and_stats(M, num_trials, sampling_freq);

        % Save per-file
        writematrix(hmm_stats, out_path);

        % Concatenate horizontally
        if isempty(all_stats)
            all_stats = hmm_stats;
        else
            % Pad shorter arrays if needed
            maxRows = max(size(all_stats,1), size(hmm_stats,1));
            % Expand all_stats
            if size(all_stats,1) < maxRows
                all_stats(end+1:maxRows, :) = NaN;
            end
            % Expand hmm_stats
            if size(hmm_stats,1) < maxRows
                hmm_stats(end+1:maxRows, :) = NaN;
            end
            % Append new block of columns
            all_stats = [all_stats, hmm_stats];
        end

    catch ME
        warning('Skipping %s due to error: %s', files(k).name, ME.message);
        continue;
    end
end

% Save concatenated results
concat_out = fullfile(output_path, 'HMMStats_all.csv');
writematrix(all_stats, concat_out);
fprintf('Saved concatenated results (horizontally stacked): %s\n', concat_out);

fprintf('Done. Processed %d file(s).\n', numel(files));
end
