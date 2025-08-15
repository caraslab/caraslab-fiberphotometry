function caraslab_preprocess_FPdata(Savedir, sel, tranges, previous_trange, select_trange, subtract_405, do_airPLS)
    % This function takes fiber photometry csv files and employs in this order:
    % 1. Low-pass filter
    % 2. Auto-detect LED onset by derivative or prompt user to select
    % ranges
    % 3. airPLS algorithm to flatten photobleaching decay
    % 4. Fit 405 onto 465 and output df/f
    % 5. Saves a filename_dff.csv file
    
    %Input variables:
    %
    %       Savedir: path to folder containing data directories. Each directory
    %                should contain a binary (-dat) data file and
    %                a kilosort configuration (config.mat) file. 
    %
    %       sel:    if 0, program will cycle through all folders
    %               in the data directory.    
    %
    %               if 1, program will prompt user to select folder
    
    %       tranges:
    %               array with the recording time range for analysis, 
    %               [0, Inf] if whole range is desired
    
    %       guess_t1:    
    %               if 1, LED onset will be determined by second derivative
    %               if 0, an attempt to read t1 will be made from ops,
    %                   if none is found, it will prompt user for selection
   
    %       select_trange:
    %               if 1, user will be prompted with a plot to select
    %               desired ranges for analysis (multiple ranges can be
    %               selected
    %               
    %       subtract_405:    
    %               if 1, 405 nm signal will be fitted then subtracted from
    %               the 465 nm signal

    %Written by M Macedo-Lima 10/05/20
    if ~sel
        datafolders = caraslab_lsdir(Savedir);
        datafolders = {datafolders.name};

    elseif sel  
        %Prompt user to select folder
        datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
        datafolders = {};
        for i=1:length(datafolders_names)
            [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
        end
    %     [~,name] = fileparts(pname);
    %     datafolders = {name};  
    %     
    end
    

    
    % Diagnostic plot colors
    global color_isosbestic color_main
    color_isosbestic = [179, 0, 179]/255;
    color_main = [0, 128, 0]/255;

    %For each data folder...
    for i = 1:numel(datafolders)
        close all;
        cur_path.name = datafolders{i};
        cur_savedir = fullfile(Savedir, cur_path.name);

        %Load in info file
        % Catch error if -mat file is not found
        try
            cur_infofiledir = dir(fullfile(cur_savedir, [cur_path.name '.info']));
            load(fullfile(cur_infofiledir.folder, cur_infofiledir.name), '-mat', 'epData');
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('\n-mat file not found\n')
                continue
            else
                fprintf(ME.identifier)
                fprintf(ME.message)
                continue
            end
        end
        fprintf('\n======================================================\n')
        fprintf('Processing and fitting photometry data, %s.......\n', cur_path.name)

        t0 = tic;

        % Create a configs variable to hold some info about the recording.
        % Optionally read a previously saved config to read info about T1 and
        % T2
        if ~previous_trange
            ops = struct();
        else
            %Load in configuration file (contains ops struct)
            % Catch error if -mat file is not found
            try
                load(fullfile(cur_savedir, 'config.mat'));
            catch ME
                if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                    fprintf('\nconfig file not found. Assuming new recording...\n')
                    ops = struct();
                else
                    fprintf(ME.identifier)
                    fprintf(ME.message)
                    continue
                end
            end        
        end

        % Load CSV file
        cur_datafiledir = dir(fullfile(cur_savedir, ['*' cur_path.name '_rawRecording.csv']));
        fullpath = fullfile(cur_datafiledir.folder, cur_datafiledir.name);
        cur_data = readtable(fullpath);

        % Use time column to calculate fs
        % Alternatively, read .info file and this should be there
        fs = 1/mean(diff(cur_data.Time));  

        ops.fraw = fullpath;
        ops.fs = fs;

        %% Moving-mean filter and remove bad recording snippets
        signal_main = cur_data.Ch465_mV;
        signal_isosbestic = cur_data.Ch405_mV;
        time_vec = cur_data.Time;
        
        % % 100-point zero-phase moving-mean filter
        % mov_mean_filter = 1/100*ones(100,1);
        % signal_main = filtfilt(mov_mean_filter, 1, signal_main);
        % signal_isosbestic = filtfilt(mov_mean_filter, 1, signal_isosbestic);
        
        % New tweakable filter
        fslow = 3;  % 3 Hz is recommended by Keevers and Bressel, 2025 but can be tweaked
        signal_main = doublefilter_FPdata(signal_main, ops.fs, fslow);
        signal_isosbestic = doublefilter_FPdata(signal_isosbestic, ops.fs, fslow);

        % Loop through the time ranges and eliminate data points outside of
        % them
        signal_main_offset = [];
        signal_isosbestic_offset = [];
        time_vec_offset = [];

        % Grab tranges from ops
        if previous_trange
            if isfield(ops, 'tranges')
                tranges = ops.tranges;
                for trange_idx=1:length(tranges)
                    cur_trange = tranges{trange_idx};
                    % Previous version of the code allowed Inf values for
                    % trange(2). Convert them here
                    if isinf(cur_trange(2))
                        cur_trange(2) = max(time_vec);
                    end
                    idx_slice = (time_vec >= cur_trange(1)) & (time_vec <= cur_trange(2));
    
                    % Fill signals with valid input
                    signal_main_offset = [signal_main_offset; signal_main(idx_slice)];
                    signal_isosbestic_offset = [signal_isosbestic_offset; signal_isosbestic(idx_slice)];
                    time_vec_offset = [time_vec_offset; time_vec(idx_slice)];

                    % Save new tranges for future use if needed
                    tranges{trange_idx} = [cur_trange(1) cur_trange(2)];
                end
            else
                fprintf('Previous time ranges not found in ops. Switching to manual selection...')
                select_trange = 1;
                previous_trange = 0;
            end
        end

        % Select from plot manually and override what's in ops
        if ~previous_trange  % Re-check here in case it's overridden 
            if select_trange
                tranges = select_tranges_from_plot(time_vec, signal_isosbestic, signal_main);
                
                for trange_idx=1:length(tranges)
                    cur_trange = tranges{trange_idx};
                    idx_slice = (time_vec >= cur_trange(1)) & (time_vec <= cur_trange(2));
    
                    % Fill signals with valid input
                    signal_main_offset = [signal_main_offset; signal_main(idx_slice)];
                    signal_isosbestic_offset = [signal_isosbestic_offset; signal_isosbestic(idx_slice)];
                    time_vec_offset = [time_vec_offset; time_vec(idx_slice)];
    
                end
                % Or detect LED onset automatically; 
                % Detect LED onset via second-derivative > 1 and eliminate everything
                % from that point + 1/fs*10 s
                % This was determined via visual inspection. Make sure to inspect!
                % Check if T1 and T2 already exist
                % NOTE: This is not set up to detect unplugging artifacts yet
            else
                for trange_idx=1:length(tranges)
                    trange = tranges{trange_idx};
    
                    % Calculate first timestamp 
                    if trange(1) == 0 && previous_trange == 1  
                        diff_465 = diff(diff(signal_main));
    
                        diff_mean = mean(diff_465);
                        diff_std = std(diff_465);
                        diff_thresh = diff_mean + 20*diff_std;
    
                        % DEBUG
        %                 figure;
        %                 plot(diff_465);
    
                        crossing = find(diff_465 > diff_thresh, 1, 'first');
                        
                        % Eliminate from start until crossing + fs*10
                        T1_idx = floor(crossing + fs*10);
                    else
                        T1_idx = max([1, floor(trange(1)*fs)]);
                    end
                    T1 = round(T1_idx/fs, 4);
                    idx_slice = (time_vec >= T1);
    
                    % Fill signals with valid input
                    signal_main_offset = [signal_main_offset; signal_main(idx_slice)];
                    signal_isosbestic_offset = [signal_isosbestic_offset; signal_isosbestic(idx_slice)];
                    time_vec_offset = [time_vec_offset; time_vec(idx_slice)];
    
                    % Save new tranges for future use if needed
                    tranges{end+1} = [T1 max(time_vec)];
                end
            end
        end

        % Update ops
        ops.tranges = tranges;
        
        %% airPLS algorithm to correct baseline
        % Baseline correction using adaptive iteratively reweighted Penalized Least Squares;		
        % These parameters seem to work well for the most part: num2cell([10e8, 1, 0.1, 0.5, 50])
        % TODO: run airPLS separately for each continuous chunk of
        % recording, i.e., stop it during large artifacts
        if do_airPLS
            airPLSconfig.input = num2cell([10e7, 1, 0.1, 0.5, 50]);
            signal_main_offset_pls = [];
            signal_main_offset_base = [];
            signal_isosbestic_offset_pls = [];
            signal_isosbestic_offset_base = [];
            for i=1:length(ops.tranges)
                idx_slice = (time_vec_offset >= ops.tranges{i}(1)) & (time_vec_offset <= ops.tranges{i}(2));

                cur_signal_main = signal_main_offset(idx_slice);
                cur_signal_isosbestic = signal_isosbestic_offset(idx_slice);

                % 465
                [temp_offset_pls, temp_offset_base]= airPLS(cur_signal_main', airPLSconfig.input{:});
                signal_main_offset_pls = [signal_main_offset_pls; temp_offset_pls'];
                signal_main_offset_base = [signal_main_offset_base; temp_offset_base'];
                
                % Isosbestic
                [temp_offset_pls, temp_offset_base]= airPLS(cur_signal_isosbestic', airPLSconfig.input{:});
                signal_isosbestic_offset_pls = [signal_isosbestic_offset_pls; temp_offset_pls'];
                signal_isosbestic_offset_base = [signal_isosbestic_offset_base; temp_offset_base'];
            end
        else
            signal_main_offset_pls = signal_main_offset;
            signal_main_offset_base = zeros(size(signal_main_offset_pls));

            signal_isosbestic_offset_pls = signal_isosbestic_offset;
            signal_isosbestic_offset_base = zeros(size(signal_isosbestic_offset_pls));
        end
        
        %% Standardize signals and regress to extract df/f
        signal_main_offset_pls = (signal_main_offset_pls - median(signal_main_offset_pls)) / std(signal_main_offset_pls);
        signal_isosbestic_offset_pls = (signal_isosbestic_offset_pls - median(signal_isosbestic_offset_pls)) / std(signal_isosbestic_offset_pls);

        %% Null-Z normalization
        % Recommended by Keevers and Bressel, 2025 for IRLS, but doesn't
        % change much
        % signal_main_offset_pls = signal_main_offset_pls ./ sqrt(mean(signal_main_offset_pls.^2));
        % signal_isosbestic_offset_pls = signal_isosbestic_offset_pls ./ sqrt(mean(signal_isosbestic_offset_pls.^2));

        %% Non-negative robust linear regression
        % If it fails, downsampling helps. Keep track of how many
        % downsampling points were required for the fit to work, then
        % reupsample everything by the same factor
%         fit_success = 0;
%         N = 1;
%         while ~fit_success && N < 100
%             [fit_result, ~, fit_success] = fit(signal_isosbestic_offset_pls, signal_main_offset_pls, ...
%                 fittype('poly1'),'Robust','off', 'lower', [0 -Inf], 'Normalize', 'off');
% 
%             fit_success = fit_success.exitflag;
% 
%             if ~fit_success
%                 signal_main_offset = downsample_sig(signal_main_offset, N);
%                 signal_isosbestic_offset = downsample_sig(signal_isosbestic_offset, N);  
%                 signal_main_offset_pls = downsample_sig(signal_main_offset_pls, N);
%                 signal_isosbestic_offset_pls = downsample_sig(signal_isosbestic_offset_pls, N);  
%                 signal_main_offset_base = downsample_sig(signal_main_offset_base, N);
%                 signal_isosbestic_offset_base = downsample_sig(signal_isosbestic_offset_base, N);  
%                 time_vec_offset = downsample_sig(time_vec_offset, N);
%                 N = N+1;
%             end
%         end
% 
%         if N == 100
%             % If fit fails after 100 attempts, you should make sure both of
%             % your signals (465 and 405) are in good quality; also check
%             % that all noise has been removed from analysis
%             fprintf('Downsampling failed after 100 attempts. Exiting...')
%             return
%         end
% 
        % Y_fit_all = fit_result(signal_isosbestic_offset_pls);
% 
%         % Simple fitting; similar to GuPPY
% %         fit_result = polyfit(signal_isosbestic_offset_pls, signal_main_offset_pls, 1);
% 
% %         Y_fit_all = (fit_result(1)*signal_isosbestic_offset_pls)+fit_result(2);
% 
%         % Upsample everything back
%         signal_main_offset = upsample_sig(signal_main_offset, N);
%         signal_isosbestic_offset = upsample_sig(signal_isosbestic_offset, N);  
%         signal_main_offset_pls = upsample_sig(signal_main_offset_pls, N);
%         signal_isosbestic_offset_pls = upsample_sig(signal_isosbestic_offset_pls, N);  
%         signal_main_offset_base = upsample_sig(signal_main_offset_base, N);
%         signal_isosbestic_offset_base = upsample_sig(signal_isosbestic_offset_base, N);
%         Y_fit_all = upsample_sig(Y_fit_all, N); 
%         time_vec_offset = upsample_sig(time_vec_offset, N);

        %% IRLS regression
        % Keevers and Bressel, 2025
        % 4.685 is the default
        % [1.4 3 4.685] were tested in the paper with slightly better
        % results
        % Tuning constants for iteratively-reweighted-least-squares regression
        % Smaller values = more aggressive downweighting of outliers 
        IRLS_constant = 1.4; 
        [~, Y_fit_all] = IRLS_dFF(signal_main_offset_pls, signal_isosbestic_offset_pls, IRLS_constant);

        %% Subtract Y_fit to get the residual 'transients'; 
        % the standardization was already performed before fitting
        if subtract_405
            signal_main_sub = signal_main_offset_pls - Y_fit_all;
        else
            signal_main_sub = signal_main_offset_pls;
        end
        %% Compile and save table
        datafilepath = split(cur_datafiledir.folder, filesep);
        subj_id = split(datafilepath{end-1}, '-');
        subj_id = join(subj_id(1:3), "-");
        datafilename = fullfile(cur_datafiledir.folder, [subj_id{1} '_' datafilepath{end}]);

        TT = array2table([time_vec_offset signal_main_offset_pls ... 
            signal_isosbestic_offset_pls Y_fit_all signal_main_sub],...
            'VariableNames',{'Time' 'Ch465_mV' 'Ch405_mV' 'Ch405_fit' 'Ch465_dff'});
        writetable(TT, [datafilename '_dff.csv']);    

        %% Plots

        close all;
        f = figure;

        % Raw recording plot
        subplot(4, 1, 1)
        plot(time_vec, signal_isosbestic, 'color', color_isosbestic, 'LineWidth', 1); hold on;
        plot(time_vec, signal_main, 'color', color_main, 'LineWidth', 1);
        
        % Highlight valid time ranges
        fill_YY = [min([signal_main; signal_isosbestic]), max([signal_main; signal_isosbestic])];
        YY = repelem(fill_YY, 1, 2);
        fill_color = [0, 0, 153]/255;
        for trange_idx=1:length(tranges)
            trange = tranges{trange_idx};
            fill_XX = [trange(1) trange(2)];
            XX = [fill_XX, fliplr(fill_XX)];
            h = fill(XX, YY, fill_color);
            % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
            set(h,'facealpha',.05,'edgecolor','none')
        end

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('mV', 'FontSize', 12)
        title(sprintf('Trial raw recording'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        % Make a legend
        legend('405 nm','465 nm','AutoUpdate', 'off');

        % airPLS baseline plots
        if do_airPLS
            subplot(8, 2, 5)
            plot(time_vec_offset, signal_main_offset, 'color', color_main, 'LineWidth', 1);  hold on;
            plot(time_vec_offset, signal_main_offset_base, 'color', 'black', 'LineWidth', 1);   
            % Finish up the plot
            axis tight
            ylabel('mV', 'FontSize', 12)
            title(sprintf('Baseline fit'))
            set(gcf, 'Position',[100, 100, 800, 500])
            set(gca,'box','off')

            subplot(8, 2, 7)
            plot(time_vec_offset, signal_isosbestic_offset, 'color', color_isosbestic, 'LineWidth', 1);  hold on;
            plot(time_vec_offset, signal_isosbestic_offset_base, 'color', 'black', 'LineWidth', 1);   
            % Finish up the plot
            axis tight
            xlabel('Time, s','FontSize',12)
            ylabel('mV', 'FontSize', 12)
        %     title(sprintf('Trial recording (405-fitted)'))
            set(gcf, 'Position',[100, 100, 800, 500])
            set(gca,'box','off')
        end
        
        subplot(8, 2, 6)
        plot(time_vec_offset, signal_main_offset_pls, 'color', color_main, 'LineWidth', 1);
        % Finish up the plot
        axis tight
        ylabel('AU', 'FontSize', 12)
        title(sprintf('Baseline corrected'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        subplot(8, 2, 8)
        plot(time_vec_offset, signal_isosbestic_offset_pls, 'color', color_isosbestic, 'LineWidth', 1);
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('AU', 'FontSize', 12)
    %     title(sprintf('Trial recording (405-fitted)'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        % 405-fit plot
        subplot(4, 1, 3)
        plot(time_vec_offset, Y_fit_all, 'color', color_isosbestic, 'LineWidth', 1);  hold on;
        plot(time_vec_offset, signal_main_offset_pls, 'color', color_main, 'LineWidth', 1);   
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('AU', 'FontSize', 12)
        title(sprintf('Trial recording (405-fitted)'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')
        % Make a legend
    %     legend('405 nm','465 nm','Onset', 'AutoUpdate', 'off');

        % 405-subtracted plot
        subplot(4, 1, 4)

        % Add AM trial events
        try
            tt0_events = epData.epocs.TTyp.onset(epData.epocs.TTyp.data == 0);
            fill_YY = [min(signal_main_sub), max(signal_main_sub)];
            YY = repelem(fill_YY, 1, 2);
            fill_color = [163, 163, 194]/255;
            hold on;
            for event_idx=1:length(tt0_events)
                fill_XX = [tt0_events(event_idx) tt0_events(event_idx)+1];
                XX = [fill_XX, fliplr(fill_XX)];
                h = fill(XX, YY, fill_color);
                % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
                set(h,'facealpha',.5,'edgecolor','none')
            end
        catch ME
            if strcmp(ME.identifier, 'MATLAB:nonExistentField')
                if isfield(epData.epocs, 'STTL')
                    tt0_events = epData.epocs.STTL.onset;
                elseif isfield(epData.epocs, 'STL_')
                    tt0_events = epData.epocs.STL_.onset;
                end
                fill_YY = [min(signal_main_sub), max(signal_main_sub)];
                YY = repelem(fill_YY, 1, 2);
                fill_color = [163, 163, 194]/255;
                hold on;
                for event_idx=1:length(tt0_events)
                    fill_XX = [tt0_events(event_idx) tt0_events(event_idx)+1];
                    XX = [fill_XX, fliplr(fill_XX)];
                    h = fill(XX, YY, fill_color);
                    % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
                    set(h,'facealpha',.5,'edgecolor','none')
                end
                % Make a legend
                legend('AM trial', 'AutoUpdate', 'off');
                
                % L trough
                tt0_events = epData.epocs.LTR_.onset;
                fill_YY = [min(signal_main_sub), max(signal_main_sub)];
                YY = repelem(fill_YY, 1, 2);
                fill_color = [184, 51, 106]/255;
                hold on;
                for event_idx=1:length(tt0_events)
                    fill_XX = [tt0_events(event_idx) tt0_events(event_idx)+1];
                    XX = [fill_XX, fliplr(fill_XX)];
                    h = fill(XX, YY, fill_color);
                    % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
                    set(h,'facealpha',.5,'edgecolor','none')
                end
                % Make a legend
                legend('L trough', 'AutoUpdate', 'off');
                
                % R trough
                tt0_events = epData.epocs.RTR_.onset;
                fill_YY = [min(signal_main_sub), max(signal_main_sub)];
                YY = repelem(fill_YY, 1, 2);
                fill_color = [66, 122, 161]/255;
                hold on;
                for event_idx=1:length(tt0_events)
                    fill_XX = [tt0_events(event_idx) tt0_events(event_idx)+1];
                    XX = [fill_XX, fliplr(fill_XX)];
                    h = fill(XX, YY, fill_color);
                    % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
                    set(h,'facealpha',.5,'edgecolor','none')
                end
                % Make a legend
                legend('R trough', 'AutoUpdate', 'off');
            end
        end


        % Plot data on top
        plot(time_vec_offset, signal_main_sub, 'color', color_main, 'LineWidth', 1); 

        % Finish up the plot
        axis tight
        xlabel('Time (s)','FontSize',12)
        ylabel('dF/F %', 'FontSize', 12)
        title(sprintf('Trial recording (405-subtracted and drift-corrected)'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        % Make a legend
        legend('AM trial', 'AutoUpdate', 'off');

        % Save .fig and .pdf
        savefig(f, [datafilename '_trialPlot.fig'])
        set(gcf, 'PaperPositionMode', 'auto', 'renderer','Painters');
        print(gcf, '-painters', '-dpdf', '-r300', [datafilename '_trialPlot'], '-fillpage')

        %Save configuration file
        configfilename  = fullfile(cur_savedir,'config.mat');
        save(configfilename,'ops')
        fprintf('Saved configuration file: %s\n', configfilename)

        tEnd = toc(t0);
        fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

    end

    function sig = downsample_sig(sig, n_points)
        sig = interp1(1:length(sig), sig, linspace(1,length(sig), length(sig)/n_points + 1));
        sig = sig';
    end

    function sig = upsample_sig(sig, n_points)
        sig = interp1(1:length(sig), sig, linspace(1,length(sig), length(sig)*n_points + 1));
        sig = sig';
    end
end

function new_tranges = select_tranges_from_plot(time_vec, signal_isosbestic, signal_main)
    global color_isosbestic color_main
    new_tranges = {};
    % Plot raw data
    temp_f = gcf;
    set(temp_f, 'Position', get(0, 'Screensize'));
    set(temp_f, 'KeyPressFcn', @myKeyPressFcn);
    hold on;
    plot(time_vec, signal_isosbestic, 'color', color_isosbestic);
    plot(time_vec, signal_main, 'color', color_main);
    fill_YY = [min([signal_main; signal_isosbestic]), max([signal_main; signal_isosbestic])];
    YY = repelem(fill_YY, 1, 2);
    

    selection_complete = 0;
    loop_idx = 1;
    while ~selection_complete
        [cur_trange, ~] = ginput(2);
        
        if length(cur_trange) == 2
            sort(cur_trange);
            if cur_trange(2) > max(time_vec)
                cur_trange(2) = max(time_vec);
            end
            fill_color = rand(1,3);
            
            fill_XX = [cur_trange(1) cur_trange(2)];
            XX = [fill_XX, fliplr(fill_XX)];
            h = fill(XX, YY, fill_color);
            % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
            set(h,'facealpha',.05,'edgecolor','none')
            new_tranges{end+1} = cur_trange';
            loop_idx = loop_idx + 1;
        else
            selection_complete = 1;
        end
        
    end
    close(temp_f)
end