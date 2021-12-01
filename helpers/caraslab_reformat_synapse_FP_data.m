function caraslab_reformat_synapse_FP_data(Tankdir,Savedir,sel)
%epData = caras_lab_reformat_FP_synapse_data(Tankdir,Savedir,sel);
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
%       
%
%   Written by ML Caras Mar 22, 2019 
%   patched by M Macedo-Lima 9/8/20

%% tweak this if needed
CH465_NAME = 'x465A';
CH405_NAME = 'x405A';

%%

%Default to cycling through all BLOCKS
if nargin < 3
    sel = 0;   
end

%Check that tank directory exists and abort if it doesn't
if ~exist(Tankdir,'dir')
    fprintf('\n Tank directory does not exist!!\n')
    return
end


%Check if save directory exists. If it doesn't, create it now.
if ~exist(Savedir,'dir')
    [success,message,messageID] = mkdir(Savedir);
    
    %Stop if directory cannot be created, and display reason why
    if ~success
        message %#ok<*NOPRT>
        messageID
        return
    end   
end

if ~sel
    %Get a list of all BLOCKS in the tank directory
    BLOCKS = caraslab_lsdir(Tankdir);
    BLOCKNAMES = {BLOCKS.name};
    

elseif sel  
    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(Tankdir,'Select data directory');
    BLOCKNAMES = {};
    for i=1:length(datafolders_names)
        [~, BLOCKNAMES{end+1}, ~] = fileparts(datafolders_names{i});
    end
end

%Check that at least one block has been selected
if isempty(BLOCKNAMES)
    fprintf('\n No BLOCKS could be found!!\n')
    return
end


%For each block
for i = 1:numel(BLOCKNAMES)
    t0 = tic;
    cur_path.name = BLOCKNAMES{i};
    cur_savedir = fullfile(Savedir, cur_path.name);

%     infofilename = fullfile(cur_savedir, [cur_path.name '_info']); % MML edit

    %Check if datafile is already saved, and if so, ask user what to do.
%     if exist(datafilename,'file')
%         reply = input('\n Reformated file exists already. Do you want to overwrite?\n Y/N:','s');
%         
%         switch reply
%             case {'n','N','no','No','NO'}
%                 continue
%         end
%     end
    
    %Convert tank data to -mat and display elapsed time
    fullpath = fullfile(Tankdir,cur_path.name);
    fprintf('\n======================================================\n')
    fprintf('Processing photometry data, %s.......\n', cur_path.name)
    
    % Try to read file. If missing, skip to the next folder
    try
        epData = TDTbin2mat(fullpath,'TYPE',{'epocs','streams'});
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            fprintf('\nFile not found\n')
            continue
        else
            fprintf(ME.identifier)
            fprintf(ME.message)
            break
        end
    end
    
    % Get some info for naming files
    subj_id = epData.info.tankpath;
    subj_id = split(subj_id, filesep);
    subj_id = subj_id{end};
    datafilename = fullfile(cur_savedir, cur_path.name);
    infofilename = fullfile(cur_savedir, [cur_path.name '.info']);
    
    % Save a -mat file with the the entire data structure and a CSV file with
    % the raw signals with a time column to simplify sampling rate 
    % Also output behavioral timestamps as csv files
    try
        mkdir(cur_savedir);
        % Write full TDT struct to file
        save([datafilename, '.mat'], 'epData', '-v7.3')
        
        % Write data streams to CSV
        % Save Stim too for possible future use
        time_vector = linspace(0, length(epData.streams.(CH465_NAME).data) / ...
            epData.streams.(CH465_NAME).fs, length(epData.streams.(CH465_NAME).data));
        
        % Audio is sampled at different sampling rate, so interpolate to
        % match data sampling
%         stim_stream_downsampled = interp1(1:length(epData.streams.Stim.data),...
%             epData.streams.Stim.data, linspace(1, length(epData.streams.Stim.data), length(time_vector)));
%         
        % Compile and save table
        fprintf('\nSaving raw streams...\n')
        
        % If audio is desired
%         TT = array2table([time_vector' epData.streams.(CH465_NAME).data' ... 
%             epData.streams.(CH405_NAME).data' stim_stream_downsampled'],...
%             'VariableNames',{'Time' 'Ch465_mV' 'Ch405_mV' 'Audio_stim'});
%         writetable(TT, [datafilename '_rawRecording.csv']);

        TT = array2table([time_vector' epData.streams.(CH465_NAME).data' ... 
            epData.streams.(CH405_NAME).data'],...
            'VariableNames',{'Time' 'Ch465_mV' 'Ch405_mV'});
        writetable(TT, [datafilename '_rawRecording.csv']);

        % Remove data streams from structure and save the rest as .info
        epData.streams.(CH465_NAME).data = [];
        epData.streams.(CH405_NAME).data = [];
%         epData.streams.Stim.data = [];
        
        fprintf('\nSaving supporting information...\n')
        save(infofilename,'epData','-v7.3')
    catch
        warning('\n ** Could not save file **.\n')
        keyboard
    end
    
    tEnd = toc(t0);
    fprintf('\n~~~~~~\nFinished in: %d minutes and %f seconds\n~~~~~~\n', floor(tEnd/60),rem(tEnd,60));


end


fprintf('\n\n ##### Finished reformatting and saving data files.\n\n')


end






