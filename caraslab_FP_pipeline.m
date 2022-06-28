%caraslab_FP_pipeline.m
%
%This pipeline transforms fiber photometry raw -sev data into -mat format, 
% processes it and outputs event timestamps and some exploratory plots
 
% Note that this pipeline was designed to be modular, i.e. you run one bit
% at a time and you can add/remove/replace modules 

% Written by M Macedo-Lima 6/21/21

%% Set your paths

% Tankdir: Where your TDT tank raw files are; This path should be a
%   subject's folder with subfolder representing different recording sessions

% Savedir: Where you wish to save the processed files

% Behaviordir: Where the ePsych behavior files are; -mat files will be
%   combined into a single file. Before running this, group similar sessions
%   into folders named: 
%   shock_training, psych_testing, pre_passive, post_passive

% chanMapSavedir: where your channel maps are

% sel: whether you want to run all or a subset of the folders. If 1, you
%   will be prompted to select folders. Multiple folders can be selected
%   using Ctrl or Shift

% rootH: path for temp Kilosort file. Should be a fast SSD


Behaviordir = '/mnt/CL_4TB_2/Matt/Fiber photometry';

% hSyn pilots
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/OFC-hSyn/SUBJ-ID-104-210608-104954';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/OFC-hSyn/SUBJ-ID-105-210608-115446';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn/SUBJ-ID-106-210608-111503';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn/SUBJ-ID-108-210616-115327';
% subtract_405 = 1;

% GRABNE pilots
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn-GRABNE/SUBJ-ID-239-211130-110221';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn-GRABNE/SUBJ-ID-240-211124-153521';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn-GRABNE/SUBJ-ID-205-211117-111337';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn-GRABNE/SUBJ-ID-206-211117-103731';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-hSyn-GRABNE/SUBJ-ID-207-211117-110621';
% subtract_405 = 0;

% AAVrg in ACx, fiber in OFC 
% Cohort 1
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-276-211223-100519';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-277';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-278-211223-104035';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-279-211223-110906';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-287-211223-115806';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-288-220113-115516';

% Cohort 2
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-336-220420-151310';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-337-220420-153301';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-338-220420-143103';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-341-220420-142655';
% Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-GCaMP8s_OFC-VO-fiber/SUBJ-ID-342-220420-140053';
Tankdir = '/mnt/CL_4TB_2/Matt/Fiber photometry/ACx-AAVrg-EGFP_OFC-VO-fiber/SUBJ-ID-343-220420-145841';
subtract_405 = 1;

Savedir =  Tankdir;

sel = 1;  % Select subfolders; 0 will run all subfolders

%% 2. CONVERT *.SEV AND TANK DATA TO *.MAT FILE
%   Function to reformat and save ephys data from TDT.
%
%   FIRST you must manually copy the .sev files from RS4 data streamer into
%   the appropriate TDT tank block.
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
%       (1) A -mat file containing an MxN matrix of raw voltages, where 
%               M = the number of channels
%               N = the number of samples
%
%       (2) A -info file containing supporting information, including
%               sampling rate, epocs, and timing
caraslab_reformat_synapse_FP_data(Tankdir,Savedir,sel);

%% 3. Output timestamps info
% This pipeline takes ePsych .mat behavioral files, combines and analyzes them and
% outputs files ready for further behavioral analyses and for aligning
% timestamps with neural recordings
% This pipeline also incorporates ephys recordings
% in the processing to extract timestamps related to spout and stimulus
% delivery

% IMPORTANT: if behavior is relevant, run this now so that createconfig can
% extract information about how much of the beginning of the recording to 
% skip due to noise

% IMPORTANT 2: organize your behavior files into subfolders to be analyzed together , e.g.
% shockTraining_pre, shockTraining_active, psychTesting_active, psychTesting_muscimol etc
% select those folders when prompted (you can select multiple folders)
caraslab_behav_pipeline(Savedir, Behaviordir, 'synapse');

%% REMOVE ARTIFACTS AND FILTER
% This function takes extracted photometry signals and processes them in
% this order:
% 1.
% 2.
% 3.
% 3

% Used to ignore points of the signal before the fitting; 
% If T1=0, automatic detection of LED onset will be performed; T2=Inf will
% read from T1 until the end
tranges = {[0 Inf]};
guess_t1 =1;
select_trange = 1;
caraslab_preprocess_FPdata(Savedir, sel, tranges, guess_t1, select_trange, subtract_405)

%% Output timestamped waves and AUCs
caraslab_getTrialData_FPdata(Savedir, sel)

%% Plot separate responses by AMdepth
only_response = {'all', 'hit', 'miss'};  % all, miss, hit

for resp_idx=1:length(only_response)
    cur_resp = only_response{resp_idx};
    caraslab_getAMDepthData_FPdata(Savedir, sel, cur_resp)
end


%% Compile data
compile_FPdata_for_analyses(Savedir)
