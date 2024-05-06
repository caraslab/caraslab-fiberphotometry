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

% Tankdir: Where your TDT tank raw files are; This path should be a
%   subject's main folder containing multiple sessions inside

% Savedir: Where you wish to save the processed files

% subtract_405: whether your 405 channel is important for removing
% movement artifacts

% sel: whether you want to run all or a subset of the folders. If 1, you
%   will be prompted to select folders. Multiple folders can be selected
%   using Ctrl or Shift

Behaviordir = '/mnt/CL_8TB_3/Matheus/VTA_FP_1IFC/TH-Cre_flex-GCaMP8m/matlab_data_files';

Tankdir = '/mnt/CL_8TB_3/temp_tank/SUBJ-ID-683-240404-133151';
Savedir =  '/mnt/CL_8TB_3/Matheus/VTA_FP_1IFC/TH-Cre_flex-GCaMP8m/SUBJ-ID-683-240404-133151';

% Tankdir = '/mnt/CL_8TB_3/temp_tank/SUBJ-ID-684-240326-150020';
% Savedir =  '/mnt/CL_8TB_3/Matheus/VTA_FP_1IFC/TH-Cre_flex-GCaMP8m/SUBJ-ID-684-240326-150020';

subtract_405 = 1;

sel = 1;  % Select subfolders; 0 will run all subfolders

%% 2. CONVERT *.SEV AND TANK DATA TO *.MAT FILE
%   Function to reformat and save data from TDT.
%
%   Input variables:
%       Tankdir:    path to tank directory
%
%       Savedir:    path to directory where -mat and -csv files will be saved
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

caraslab_reformat_synapse_FP_data(Tankdir,Savedir,sel);

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
caraslab_behav_pipeline(Savedir, Behaviordir, 'synapse_1IFC');

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
guess_t1 =1;  % If guess_t1=1, previous T1 identification saved in config will not be assumed
select_trange = 0;  % Prompt user for selecting time range for analysis
caraslab_preprocess_FPdata(Savedir, sel, tranges, guess_t1, select_trange, subtract_405)

%%
caraslab_get1IFCTrialData_FPdata(Savedir, sel)

%% Code below is specific to Caras Lab aversive AM detection paradigm
%% 5. Output timestamped waves and AUCs
caraslab_getTrialData_FPdata(Savedir, sel)

%% 6. Plot separate responses by AMdepth
only_response = {'all', 'hit', 'miss'};  % all, miss, hit

for resp_idx=1:length(only_response)
    cur_resp = only_response{resp_idx};
    caraslab_getAMDepthData_FPdata(Savedir, sel, cur_resp)
end


%% Compile data
compile_FPdata_for_analyses(Savedir)
