clear all;close all;clc

% Define the file to process
filepath = 'example_data/young_treated_trial_1.mat';

% Define parameters
ch_height  = 22;          % Channel height
De_np      = 27.65;         % Calibrated De_np of the device
wC         = 10.5;          % Contraction pore width
thresholds = [2e-5,2e-4];   % Pulse detection threshold
sampleRate = 50000;         % DAQ sampling frequency

% Call the processing function
output_table = mNPS_procData(filepath, ch_height, De_np, wC, thresholds, sampleRate);
%%
% Output
[savePath, saveName] = fileparts(filepath);
[file, path] = uiputfile('*.xlsx', 'Save output table as', ...
    fullfile(savePath, [saveName '_output.xlsx']));
if isequal(file, 0)
    fprintf('Save cancelled.\n');
    return
end
excelFile = fullfile(path, file);

writetable(output_table, excelFile);
fprintf('Saved %d cell(s) to %s\n', height(output_table), excelFile);
