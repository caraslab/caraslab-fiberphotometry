function compile_FPdata_for_analyses(Savedir)
%
% This function loops through recording folders and extracts all relevant
% files into a folder called Data inside the parent directory. The purpose
% is to centralize all subjects' data into common directories
%


%Prompt user to select folder
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end


% Create data analysis folders in parent directory
full_path_split = split(Savedir, filesep);
parent_path = strjoin(full_path_split(1:end-1), filesep);
mkdir(fullfile(parent_path, 'Data'));

behavioralPerformance_path = fullfile(parent_path, 'Data', 'Behavioral performance');
mkdir(behavioralPerformance_path);

keyFiles_path = fullfile(parent_path, 'Data', 'Key files');
mkdir(keyFiles_path);

% amTrialdPrime_path = fullfile(parent_path, 'Data', 'Output', 'AMTrial neural d-prime');
% mkdir(amTrialdPrime_path);

wholeSignal_path = fullfile(parent_path, 'Data', 'Whole session signal');
mkdir(wholeSignal_path);

%For each data folder...
for i = 1:numel(datafolders)
        cur_path.name = datafolders{i};
        cur_savedir = fullfile(Savedir, cur_path.name);
        
        % Copy Key files
        filedirs = dir(fullfile(cur_savedir, 'Info files', '*trialInfo.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(keyFiles_path, cur_filedir.name));
            end
        end
        filedirs = dir(fullfile(cur_savedir, 'Info files', '*spoutTimestamps.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(keyFiles_path, cur_filedir.name));
            end
        end
        
        % % Copy AMTrial neural d-prime
        % filedirs = dir(fullfile(cur_savedir, '*_neuralDprime.csv'));
        % for file_idx=1:length(filedirs)
        %     if ~isempty(filedirs)
        %         cur_filedir = filedirs(file_idx);
        %         copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(amTrialdPrime_path, cur_filedir.name));
        %     end
        % end 
        
        % Copy Whole session dff z-score
        filedirs = dir(fullfile(cur_savedir, '*_dff.csv'));
        for file_idx=1:length(filedirs)
            if ~isempty(filedirs)
                cur_filedir = filedirs(file_idx);
                copyfile(fullfile(cur_filedir.folder, cur_filedir.name), fullfile(wholeSignal_path, cur_filedir.name));
            end
        end
end

% Lastly copy behavior files
filedirs = caraslab_lsdir(fullfile(Savedir, 'Behavior'));
filedirs = {filedirs.name};
for file_idx=1:length(filedirs)
    if ~isempty(filedirs)
        cur_filedir = fullfile(Savedir, 'Behavior', filedirs{file_idx});
        copyfile(cur_filedir, fullfile(behavioralPerformance_path, filedirs{file_idx}));
    end
end