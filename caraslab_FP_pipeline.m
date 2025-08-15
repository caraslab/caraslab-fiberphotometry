%caraslab_FP_pipeline.m
%
%This pipeline transforms fiber photometry raw -sev data into -mat format, 
% processes it and outputs event timestamps and some exploratory plots
 
% Note that this pipeline was designed to be modular, i.e. you run one bit
% at a time and you can add/remove/replace modules 

% Written by M Macedo-Lima 6/21/21

%% 1. Set your paths
% Behaviordir: Where the ePsych behavior files are; -mat files will be
%   combined into a single file. Before running this, group similar sessions
%   into folders named: 
%   shock_training, psych_testing, pre_passive, post_passive

% input_dir: Where your TDT tank raw files are; This path should be a
%   subject's main folder containing multiple sessions inside

% output_dir: Where you wish to save the processed files

% subtract_405: whether your 405 channel is important for removing
% movement artifacts

% sel: whether you want to run all or a subset of the folders. If 1, you
%   will be prompted to select folders. Multiple folders can be selected
%   using Ctrl or Shift

% % % 
clear;

% Set paths
behavior_dir = '/path1/matlab_data_files';
root_tank_dir = '/path2';  % Input 
output_dir = '/path3';  % Output

% Uncomment one at a time
% tank_id = 'SUBJ-ID-718-240605-094128';
% tank_id = 'SUBJ-ID-719-240605-100955';
% tank_id = 'SUBJ-ID-938-250428-175742';
tank_id = 'SUBJ-ID-1011-250710-164124';

input_dir = fullfile(root_tank_dir, tank_id);
output_dir =  fullfile(output_dir, tank_id); 

subtract_405 = 1;

sel = 1;  % Select subfolders; 0 will run all subfolders

%% 2. CONVERT *.SEV AND TANK DATA TO *.MAT FILE
%   Function to reformat and save data from TDT.
%
%   Input variables:
%       input_dir:    path to tank directory
%
%       output_dir:    path to directory where -mat and -csv files will be saved
%
%       sel:        if 0 or omitted, program will cycle through all BLOCKS
%                   in the tank directory. 
%                   
%                   if 1, user will be prompted to select a BLOCK
%  
%   Uses TDTbin2mat to reformat tank data to a matlab struct. 
%   Two files are saved:    
%       (1) A -csv file containing the raw photometry data, including 465
%       and 405 nm channels and a Time column
%
%       (2) A -info file containing supporting information, including
%               sampling rate and epocs
CH465_NAME = 'x465A';
CH405_NAME = 'x405A';  % Rig 1
% CH405_NAME = 'x415A';  % Rig 3
caraslab_reformat_synapse_FP_data(input_dir,output_dir,sel, CH465_NAME, CH405_NAME);

%% 3. Output timestamps info
% This pipeline takes ePsych .mat behavioral files, combines and analyzes them and
% outputs files ready for further behavioral analyses and for aligning
% timestamps with neural recordings
% This pipeline also incorporates ephys recordings
% in the processing to extract timestamps related to spout and stimulus
% delivery

% IMPORTANT: organize your behavior files into subfolders to be analyzed together , e.g.
% shockTraining_pre, shockTraining_active, psychTesting_active, psychTesting_muscimol etc
% select those folders when prompted (you can select multiple folders)

% experiment_type: optional: 'synapse', 'intan', 'optoBehavior', '1IFC', 'synapse_1IFC'
% All types without the 1IFC tag assume aversive AM detection
experiment_type = 'synapse';
caraslab_behav_pipeline(output_dir, behavior_dir, 'experiment_type', experiment_type);


%% 4. REMOVE ARTIFACTS AND FILTER
% This function takes extracted photometry signals and processes them in
% this order:
% 1. Low-pass filter
% 2. Auto-detect LED onset by derivative or prompt user to select
% ranges
% 3. airPLS algorithm to flatten photobleaching decay
% 4. Fit 405 onto 465 and output df/f
% 5. Saves a filename_dff.csv file

% Ignore points of the signal before the fitting; 
% If T1=0, automatic detection of LED onset will be performed; T2=Inf will
% read from T1 until the end
tranges = {[0 Inf]};
previous_trange = 1;  % If guess_t1=1, previous T1 identification saved in config will not be assumed
select_trange = 1;  % Prompt user for selecting time range for analysis
do_airPLS = 1;
caraslab_preprocess_FPdata(output_dir, sel, tranges, previous_trange, select_trange, subtract_405, do_airPLS)

%% Steps 5-7 are optional and for quick visualization only and will not impact analysis

%% 5. Output timestamped waves and AUCs for 1IFC protocol
caraslab_get1IFCTrialData_FPdata(output_dir, sel)

%% 6. Output timestamped waves and AUCs for AversiveAM protocol
caraslab_getTrialData_FPdata(output_dir, sel)

%% 7. Plot separate responses by AMdepth
only_response = {'all'};  % all, miss, hit

for resp_idx=1:length(only_response)
    cur_resp = only_response{resp_idx};
    caraslab_getAMDepthData_FPdata(output_dir, sel, cur_resp)
end


%% 8. Compile data for Python pipeline
compile_FPdata_for_analyses(output_dir)
