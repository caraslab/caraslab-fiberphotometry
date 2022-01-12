function caraslab_preprocess_FPdata(Savedir, sel, tranges, guess_t1)
    % This function takes fiber photometry csv files and employs in this order:
    % 1. Downsampling 10x by interpolation
    % 2. Low-pass filter
    % 3. Upsampling 10x by interpolation
    % 4. Auto-detect LED onset by derivative
    % 5. Fit 405 onto 465 and output df/f
    % 4. Saves a filename_dff.csv file
    %
    %Input variables:
    %
    %       Savedir: path to folder containing data directories. Each directory
    %                should contain a binary (-dat) data file and
    %                a kilosort configuration (config.mat) file. 
    %
    %       sel:    if 0 or omitted, program will cycle through all folders
    %               in the data directory.    
    %
    %               if 1, program will prompt user to select folder

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


    %For each data folder...
    for i = 1:numel(datafolders)

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

        % Moving-mean, filter and normalize by isosbestic signal
        % Also interpolate time vector to keep timepoint numbers consistent
        % Add audio stream here too when needed
        y_465 = cur_data.Ch465_mV;
        y_405 = cur_data.Ch405_mV;
        time_vec = cur_data.Time;
        
        % Filter
        % Downsample here so fitting improves and everything speeds  up

        % Pad with "0"s to avoid filter artifacts
%         padsize = 1000;
%         y_465 = [repmat(y_465(1,:), [padsize 1]); y_465; repmat(y_465(end,:), [padsize 1])];
%         y_405 = [repmat(y_405(1,:), [padsize 1]); y_405; repmat(y_405(end,:), [padsize 1])];

%         lowpass = 5;
% %         [b1, a1] = butter(3, 2*lowpass/(fs/N), 'low'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
%         [b1, a1] = butter(3, 2*lowpass/(fs), 'low'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
% 
%         ops.lowpass = lowpass;
% 
%         % Filter then remove padding
%         y_465 = filter(b1, a1, y_465); % causal forward filter
%         y_465 = flipud(y_465); % reverse time
%         y_465 = filter(b1, a1, y_465); % causal forward filter again
%         y_465 = flipud(y_465); % reverse time back
% 
%         y_405 = filter(b1, a1, y_405); % causal forward filter
%         y_405 = flipud(y_405); % reverse time
%         y_405 = filter(b1, a1, y_405); % causal forward filter again
%         y_405 = flipud(y_405); % reverse time back
% 
%         y_465 = y_465(padsize+1:end-padsize, :);
%         y_405 = y_405(padsize+1:end-padsize, :);

        % Zero-phase moving-mean
        mov_mean_filter = 1/100*ones(100,1);
        y_465 = filtfilt(mov_mean_filter, 1, y_465);
        y_405 = filtfilt(mov_mean_filter, 1, y_405);

        % Detect LED onset via first-derivative > 1 and eliminate everything
        % from that point + 30 s
        % This was determined via visual inspection. Make sure to inspect!
        % Check if T1 and T2 already exist
        if ~guess_t1
            if isfield(ops, 'tranges')
                tranges = ops.tranges;
            end
        end
        
        % Loop through the time ranges and eliminate data points outside of
        % them
        y_465_offset = [];
        y_405_offset = [];
        time_vec_offset = [];
        new_tranges = {};
        for trange_idx=1:length(tranges)
            trange = tranges{trange_idx};
            
            % Calculate first timestamp 
            if trange(1) == 0  && guess_t1 == 1  
                diff_thresh = 1.5;
                diff_465 = diff(y_465);
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
                T2_idx = min(length(y_465), ceil(trange(2)*fs));
            end
            
            % Fill signals with valid input
            y_465_offset = [y_465_offset; y_465(T1_idx:T2_idx)];
            y_405_offset = [y_405_offset; y_405(T1_idx:T2_idx)];
            time_vec_offset = [time_vec_offset; time_vec(T1_idx:T2_idx)];
            
            % Save new tranges for future use if needed
            new_tranges{end+1} = [T1_idx/fs T2_idx/fs];
        end
        
        % Update ops
        tranges = new_tranges;
        ops.tranges = tranges;
        
        % airPLS algorithm to correct baseline
        airPLSconfig.input = num2cell([10e9, 1, 0.1, 0.5, 50]);
        [y_465_offset_pls,y_465_offset_base]= airPLS(y_465_offset', airPLSconfig.input{:});
        y_465_offset_pls = y_465_offset_pls';
        y_465_offset_base = y_465_offset_base';
    %     y_465_offset = y_465_offset - Z';

        [y_405_offset_pls,y_405_offset_base]= airPLS(y_405_offset', airPLSconfig.input{:});
        y_405_offset_pls = y_405_offset_pls';
        y_405_offset_base = y_405_offset_base';
    %     y_405_offset = y_405_offset - Z';


        % Standardize signals
        y_465_offset_pls = (y_465_offset_pls - median(y_465_offset_pls)) / std(y_465_offset_pls);
        y_405_offset_pls = (y_405_offset_pls - median(y_405_offset_pls)) / std(y_405_offset_pls);

        % Simple regression
        % regress FP signal against 405 control to learn coeffs:
        % (if numerical problems, try subsampling data)
    %     bls = polyfit(y_405_offset_pls,y_465_offset_pls,1);
    %     Y_fit_all = bls(1) .* y_405_offset_pls + bls(2);

        % Non-negative robust linear regression; sometimes this fails; not sure
        % what to do in this case... downsampling?
        
        fit_success = 0;
        N = 1;
        while ~fit_success && N < 100
            
            % Fit function cannot handle nans so mask them out here
            [bls, ~, O] = fit(y_405_offset_pls, y_465_offset_pls, ...
                fittype('poly1'),'Robust','on', 'lower', [0 -Inf], 'Normalize', 'on');
            fit_success = O.exitflag;

            if ~fit_success
                y_465_offset = downsample_sig(y_465_offset, N);
                y_405_offset = downsample_sig(y_405_offset, N);  
                y_465_offset_pls = downsample_sig(y_465_offset_pls, N);
                y_405_offset_pls = downsample_sig(y_405_offset_pls, N);  
                y_465_offset_base = downsample_sig(y_465_offset_base, N);
                y_405_offset_base = downsample_sig(y_405_offset_base, N);  
                time_vec_offset = downsample_sig(time_vec_offset, N);
                N = N+1;
            end
        end
        fprintf('Downsampling factor needed for fitting: %d\n', N)
        
        Y_fit_all = bls(y_405_offset_pls);
        
        % Upsample everything
        y_465_offset = upsample_sig(y_465_offset, N);
        y_405_offset = upsample_sig(y_405_offset, N);  
        y_465_offset_pls = upsample_sig(y_465_offset_pls, N);
        y_405_offset_pls = upsample_sig(y_405_offset_pls, N);  
        y_465_offset_base = upsample_sig(y_465_offset_base, N);
        y_405_offset_base = upsample_sig(y_405_offset_base, N);
        Y_fit_all = upsample_sig(Y_fit_all, N); 
        time_vec_offset = upsample_sig(time_vec_offset, N);
        
        % Subtract Y_fit to get the residual 'transients'; 
        % the standardization was already performed before fitting
        y_465_sub = y_465_offset_pls - Y_fit_all;

        % Compile and save table
        datafilepath = split(cur_datafiledir.folder, filesep);
        subj_id = split(datafilepath{end-1}, '-');
        subj_id = join(subj_id(1:3), "-");
        datafilename = fullfile(cur_datafiledir.folder, [subj_id{1} '_' datafilepath{end}]);

        TT = array2table([time_vec_offset y_465_offset_pls ... 
            y_405_offset_pls Y_fit_all y_465_sub],...
            'VariableNames',{'Time' 'Ch465_mV' 'Ch405_mV' 'Ch405_fit' 'Ch465_dff'});
        writetable(TT, [datafilename '_dff.csv']);    

        % Diagnostic Plots 
        color_405 = [179, 0, 179]/255;
        color_465 = [0, 128, 0]/255;

        close all;
        f = figure;

        % Raw recording plot
        subplot(4, 1, 1)
        plot(time_vec, y_405, 'color', color_405, 'LineWidth', 1); hold on;
        plot(time_vec, y_465, 'color', color_465, 'LineWidth', 1);
        
        % Highlight valid time ranges
        fill_YY = [min(y_465), max(y_465)];
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
        legend('405 nm','465 nm','T1', 'T2', 'AutoUpdate', 'off');


        % airPLS baseline plots
        subplot(8, 2, 5)
        plot(time_vec_offset, y_465_offset, 'color', color_465, 'LineWidth', 1);  hold on;
        plot(time_vec_offset, y_465_offset_base, 'color', 'black', 'LineWidth', 1);   
        % Finish up the plot
        axis tight
        ylabel('mV', 'FontSize', 12)
        title(sprintf('Baseline fit'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        subplot(8, 2, 7)
        plot(time_vec_offset, y_405_offset, 'color', color_405, 'LineWidth', 1);  hold on;
        plot(time_vec_offset, y_405_offset_base, 'color', 'black', 'LineWidth', 1);   
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('mV', 'FontSize', 12)
    %     title(sprintf('Trial recording (405-fitted)'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        subplot(8, 2, 6)
        plot(time_vec_offset, y_465_offset_pls, 'color', color_465, 'LineWidth', 1);
        % Finish up the plot
        axis tight
        ylabel('AU', 'FontSize', 12)
        title(sprintf('Baseline corrected'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')

        subplot(8, 2, 8)
        plot(time_vec_offset, y_405_offset_pls, 'color', color_405, 'LineWidth', 1);
        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('AU', 'FontSize', 12)
    %     title(sprintf('Trial recording (405-fitted)'))
        set(gcf, 'Position',[100, 100, 800, 500])
        set(gca,'box','off')


        % 405-fit plot
        subplot(4, 1, 3)
        plot(time_vec_offset, Y_fit_all, 'color', color_405, 'LineWidth', 1);  hold on;
        plot(time_vec_offset, y_465_offset_pls, 'color', color_465, 'LineWidth', 1);   
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
        tt0_events = epData.epocs.TTyp.onset(epData.epocs.TTyp.data == 0);
        fill_YY = [min(y_465_sub), max(y_465_sub)];
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

        % Plot data on top
        plot(time_vec_offset, y_465_sub, 'color', color_465, 'LineWidth', 1); 

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
        print(gcf, '-painters', '-dpdf', '-r300', [datafilename '_recordingPlot'], '-fillpage')

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