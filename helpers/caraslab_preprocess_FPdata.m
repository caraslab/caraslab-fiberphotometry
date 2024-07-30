function caraslab_preprocess_FPdata(Savedir, sel, tranges, guess_t1, select_trange, subtract_405)
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
        if guess_t1
            ops = struct();
        else
            %Load in configuration file (contains ops struct)
            % Catch error if -mat file is not found
            try
                load(fullfile(cur_savedir, 'config.mat'));
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
        
        % 100-point zero-phase moving-mean filter
        mov_mean_filter = 1/100*ones(100,1);
        signal_main = filtfilt(mov_mean_filter, 1, signal_main);
        signal_isosbestic = filtfilt(mov_mean_filter, 1, signal_isosbestic);

        % Detect LED onset via second-derivative > 1 and eliminate everything
        % from that point + 1/fs*10 s
        % This was determined via visual inspection. Make sure to inspect!
        % Check if T1 and T2 already exist
        if ~guess_t1
            if isfield(ops, 'tranges')
                tranges = ops.tranges;
            end
        end
        
        % Loop through the time ranges and eliminate data points outside of
        % them
        signal_main_offset = [];
        signal_isosbestic_offset = [];
        time_vec_offset = [];
        new_tranges = {};
        
        if select_trange
            new_tranges = select_tranges_from_plot(time_vec, signal_isosbestic, signal_main);
            
            for trange_idx=1:length(new_tranges)
                cur_trange = new_tranges{trange_idx};
                T1_idx = round(cur_trange(1)*fs, 0);
                T2_idx = round(cur_trange(2)*fs, 0);
                % Fill signals with valid input
                if ~isinf(T2_idx)
                    signal_main_offset = [signal_main_offset; signal_main(T1_idx:T2_idx)];
                    signal_isosbestic_offset = [signal_isosbestic_offset; signal_isosbestic(T1_idx:T2_idx)];
                    time_vec_offset = [time_vec_offset; time_vec(T1_idx:T2_idx)];
                else
                    signal_main_offset = [signal_main_offset; signal_main(T1_idx:end)];
                    signal_isosbestic_offset = [signal_isosbestic_offset; signal_isosbestic(T1_idx:end)];
                    time_vec_offset = [time_vec_offset; time_vec(T1_idx:end)];
                end
            end
        else
            for trange_idx=1:length(tranges)
                trange = tranges{trange_idx};

                % Calculate first timestamp 
                if trange(1) == 0 && guess_t1 == 1  
                    diff_465 = diff(diff(signal_main));

                    diff_mean = mean(diff_465);
                    diff_std = std(diff_465);
                    diff_thresh = diff_mean + 20*diff_std;

                    % DEBUG
    %                 figure;
    %                 plot(diff_465);

                    crossing = find(diff_465 > diff_thresh, 1, 'first');
                    
                    % Eliminate from start until crossing + fs*10
                    T1_idx = round(crossing + fs*10);
                else
                    T1_idx = max([1, floor(trange(1)*fs)]);
                end

                % Set T2 to whole recording if T2==Inf
                if trange(2) == Inf
                    T2_idx = length(time_vec);
                else
                    T2_idx = min(length(signal_main), ceil(trange(2)*fs));
                end

                % Fill signals with valid input
                signal_main_offset = [signal_main_offset; signal_main(T1_idx:T2_idx)];
                signal_isosbestic_offset = [signal_isosbestic_offset; signal_isosbestic(T1_idx:T2_idx)];
                time_vec_offset = [time_vec_offset; time_vec(T1_idx:T2_idx)];

                % Save new tranges for future use if needed
                new_tranges{end+1} = [T1_idx/fs T2_idx/fs];
            end
        end

        % Update ops
        ops.tranges = new_tranges;
        
        %% airPLS algorithm to correct baseline
        % Baseline correction using adaptive iteratively reweighted Penalized Least Squares;		
        % Default parameters seem to work well for the most part: num2cell([10e8, 1, 0.1, 0.5, 50])
        airPLSconfig.input = num2cell([10e8, 1, 0.1, 0.5, 50]);
        
        [signal_main_offset_pls,signal_main_offset_base]= airPLS(signal_main_offset', airPLSconfig.input{:});
        signal_main_offset_pls = signal_main_offset_pls';
        signal_main_offset_base = signal_main_offset_base';

        [signal_isosbestic_offset_pls,signal_isosbestic_offset_base]= airPLS(signal_isosbestic_offset', airPLSconfig.input{:});
        signal_isosbestic_offset_pls = signal_isosbestic_offset_pls';
        signal_isosbestic_offset_base = signal_isosbestic_offset_base';

        %% Standardize signals and regress to extract df/f
        signal_main_offset_pls = (signal_main_offset_pls - median(signal_main_offset_pls)) / std(signal_main_offset_pls);
        signal_isosbestic_offset_pls = (signal_isosbestic_offset_pls - median(signal_isosbestic_offset_pls)) / std(signal_isosbestic_offset_pls);

        % Non-negative robust linear regression
        % If it fails, downsampling helps. Keep track of how many
        % downsampling points were required for the fit to work, then
        % reupsample everything by the same factor
        fit_success = 0;
        N = 1;
        while ~fit_success && N < 100
            [fit_result, ~, fit_success] = fit(signal_isosbestic_offset_pls, signal_main_offset_pls, ...
                fittype('poly1'),'Robust','on', 'lower', [0 -Inf], 'Normalize', 'off');
            
            fit_success = fit_success.exitflag;

            if ~fit_success
                signal_main_offset = downsample_sig(signal_main_offset, N);
                signal_isosbestic_offset = downsample_sig(signal_isosbestic_offset, N);  
                signal_main_offset_pls = downsample_sig(signal_main_offset_pls, N);
                signal_isosbestic_offset_pls = downsample_sig(signal_isosbestic_offset_pls, N);  
                signal_main_offset_base = downsample_sig(signal_main_offset_base, N);
                signal_isosbestic_offset_base = downsample_sig(signal_isosbestic_offset_base, N);  
                time_vec_offset = downsample_sig(time_vec_offset, N);
                N = N+1;
            end
        end
        
        if N == 100
            % If fit fails after 100 attempts, you should make sure both of
            % your signals (465 and 405) are in good quality; also check
            % that all noise has been removed from analysis
            fprintf('Downsampling failed after 100 attempts. Exiting...')
            return
        end
        
        Y_fit_all = fit_result(signal_isosbestic_offset_pls);

        % Simple fitting; similar to GuPPY
%         fit_result = polyfit(signal_isosbestic_offset_pls, signal_main_offset_pls, 1);
        
%         Y_fit_all = (fit_result(1)*signal_isosbestic_offset_pls)+fit_result(2);
        
        % Upsample everything back
        signal_main_offset = upsample_sig(signal_main_offset, N);
        signal_isosbestic_offset = upsample_sig(signal_isosbestic_offset, N);  
        signal_main_offset_pls = upsample_sig(signal_main_offset_pls, N);
        signal_isosbestic_offset_pls = upsample_sig(signal_isosbestic_offset_pls, N);  
        signal_main_offset_base = upsample_sig(signal_main_offset_base, N);
        signal_isosbestic_offset_base = upsample_sig(signal_isosbestic_offset_base, N);
        Y_fit_all = upsample_sig(Y_fit_all, N); 
        time_vec_offset = upsample_sig(time_vec_offset, N);
        
        % Subtract Y_fit to get the residual 'transients'; 
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
        for trange_idx=1:length(new_tranges)
            trange = new_tranges{trange_idx};
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
                cur_trange(2) = Inf;
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