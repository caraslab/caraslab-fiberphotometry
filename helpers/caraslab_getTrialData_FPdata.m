function caraslab_getTrialData_FPdata(Savedir, sel)

%% Timestamps of interest
% Change this if needed; It will skip non-present ones
% Organized in {NAME, CODE, ONSET/OFFSET}
epocs_names = {{'Spou', 1, 'onset'}, {'Spou', 1, 'offset'}, {'TTyp', 0, 'onset'}};

%% Some plotting parameters
plot_window = [-2, 6];
baseline_window = [-2, -1];

%%

if ~sel
    datafolders = caraslab_lsdir(Savedir);
    datafolders = {datafolders.name};

elseif sel  
    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
    datafolders = {};
    for dummy_idx=1:length(datafolders_names)
        [~, datafolders{end+1}, ~] = fileparts(datafolders_names{dummy_idx});
    end
%     [~,name] = fileparts(pname);
%     datafolders = {name};  
%     
end


%For each data folder...
for dummy_idx = 1:numel(datafolders)
    
    cur_path.name = datafolders{dummy_idx};
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
    fprintf('Z-scoring photometry data, %s.......\n', cur_path.name)
    
    t0 = tic;
    % Load CSV file
    cur_datafiledir = dir(fullfile(cur_savedir, ['*' cur_path.name '_dff.csv']));
    fullpath = fullfile(cur_datafiledir.folder, cur_datafiledir.name);
    cur_data = readtable(fullpath);
    time_vec = cur_data.Time;
    fs = 1/mean(diff(time_vec));
    
    % For each epoc
    for epoc_info=epocs_names
        cur_epoc_name = epoc_info{1}{1};
        cur_epoc_code = epoc_info{1}{2};
        cur_epoc_onsetoffset = epoc_info{1}{3};

        try
            code_filter = epData.epocs.(cur_epoc_name).data == cur_epoc_code;
            cur_timestamps = epData.epocs.(cur_epoc_name).(cur_epoc_onsetoffset)(code_filter);
        catch ME %#ok<NASGU>
            continue
        end
        
        % Also skip if timestamps are too few
        if length(cur_timestamps) < 3
            continue
        end
        
        time_filters = [cur_timestamps + plot_window(1) cur_timestamps + plot_window(2)];
        
        % Use this to eliminate rasters with too few points, e.g. the last
        % trial cut off by exiting ePsych
        estimated_n_points = round(diff(plot_window)*fs);
        
        raster_405 = {};
        raster_465 = {};
        raster_405fit = {};
        raster_dff = {};
        for time_slice_row=1:size(time_filters, 1)
            time_slice = time_filters(time_slice_row,:);
            index_filter = time_vec >= time_slice(1) & time_vec <= time_slice(2);
            
            % Rasters are expected to be within 1 point of each other
            % Exclude "incomplete" rasters
            if sum(index_filter) < estimated_n_points - 5
                raster_405{end+1,1} = nan(1, estimated_n_points);
                raster_465{end+1,1} = nan(1, estimated_n_points);
                raster_405fit{end+1,1} = nan(1, estimated_n_points);
                raster_dff{end+1,1} = nan(1, estimated_n_points);
            else
                raster_405{end+1,1} = cur_data.Ch405_mV(index_filter)';
                raster_465{end+1,1} = cur_data.Ch465_mV(index_filter)';
                raster_405fit{end+1,1} = cur_data.Ch405_fit(index_filter)';
                raster_dff{end+1,1} = cur_data.Ch465_dff(index_filter)';
            end
        end
        
        % Eliminate epocs with less than 5 repetitions; this helps get rid
        % of "spout" events (manual triggers) during passive
        if size(raster_465, 1) < 5
            continue
        end
        
        % Applying a time filter to a uniformly sampled signal means that the
        % length of each segment could vary by one sample.  Let's find the minimum
        % length so we can trim the excess off before calculating the mean.        
        minLength = min(cellfun('prodofsize', raster_dff));
        
        % Remove extra  point and convert to matrix
        raster_405 = cell2mat(cellfun(@(x) x(1:minLength), raster_405, 'UniformOutput',false));
        raster_465 = cell2mat(cellfun(@(x) x(1:minLength), raster_465, 'UniformOutput',false));
        raster_405fit = cell2mat(cellfun(@(x) x(1:minLength), raster_405fit, 'UniformOutput',false));
        raster_dff = cell2mat(cellfun(@(x) x(1:minLength), raster_dff, 'UniformOutput',false));

        %% %% Plot Epoch Averaged Response %% %%
        color_405 = [179/255, 0, 179/255];
        color_465 = [0, 128/255, 0];
        x_all = plot_window(1) + (1:minLength) / fs;

        %% Plot the 405 and 465 average signals dcOffsetted
        dcOff_405 = raster_405 - mean(raster_405, 'all', 'omitnan');
        dcOff_465 = raster_465 - mean(raster_465, 'all', 'omitnan');
        mean_405 = mean(dcOff_405, 'omitnan');
        mean_465 = mean(dcOff_465, 'omitnan');
        ste_405 = std(double(dcOff_405), 'omitnan')/sqrt(size(dcOff_405(all(~isnan(dcOff_405),2),:),1));
        ste_465 = std(double(dcOff_465), 'omitnan')/sqrt(size(dcOff_465(all(~isnan(dcOff_465),2),:),1));
        
        close all;
        f1 = figure;
        subplot(3,3,1)
        plot(x_all, mean_405, 'color', color_405, 'LineWidth', 3); hold on;
        plot(x_all, mean_465, 'color', color_465, 'LineWidth', 3); 

        % Make a legend
        legend('405 nm','465 nm', 'AutoUpdate', 'off');

        % Create the standard error bands for the 405 signal
        XX = [x_all, fliplr(x_all)];
        YY = [mean_405 + ste_405, fliplr(mean_405 - ste_405)];

        % Plot filled standard error bands.
        h = fill(XX, YY, color_405);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Repeat for 465
        XX = [x_all, fliplr(x_all)];
        YY = [mean_465 + ste_465, fliplr(mean_465 - ste_465)];
        h = fill(XX, YY, color_465);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('mV', 'FontSize', 12)
        title('Response (Unfitted; DC-offseted)')
        set(gcf, 'Position',[100, 100, 800, 500])

        %% Plot 465 with 405 fit
        % Get mean and ste of 405 fit signals
        mean_405 = mean(raster_405fit, 'omitnan');
        mean_465 = mean(raster_465, 'omitnan');
        ste_405 = std(double(raster_405fit), 'omitnan')/sqrt(size(raster_405fit(all(~isnan(raster_405fit),2),:),1));
        ste_465 = std(double(dcOff_465), 'omitnan')/sqrt(size(dcOff_465(all(~isnan(dcOff_465),2),:),1));
        
        subplot(3,3,2)
        plot(x_all, mean_405, 'color', color_405, 'LineWidth', 3); hold on;
        plot(x_all, mean_465, 'color', color_465, 'LineWidth', 3); 

        % Make a legend
        legend('405 nm','465 nm', 'AutoUpdate', 'off');

        % Create the standard error bands for the 405 signal
        XX = [x_all, fliplr(x_all)];
        YY = [mean_405 + ste_405, fliplr(mean_405 - ste_405)];

        % Plot filled standard error bands.
        h = fill(XX, YY, color_405);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Repeat for 465
        XX = [x_all, fliplr(x_all)];
        YY = [mean_465 + ste_465, fliplr(mean_465 - ste_465)];
        h = fill(XX, YY, color_465);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('mV', 'FontSize', 12)
        title('Response (Fitted)')
        set(gcf, 'Position',[100, 100, 800, 500])

        %% dF/F
        % Plot heat map
        subplot(3,3,4)

        imagesc(x_all, 1, raster_dff);
        colormap('jet'); % c1 = colorbar; 
        title('dF/F Heat Map (405-subtracted)');
        ylabel('Trials', 'FontSize', 12);

        % Fill band values for second subplot. Doing here to scale onset bar
        % correctly
        dff_mean = mean(raster_dff, 'omitnan');
        dff_ste = std(raster_dff, 'omitnan')/sqrt(size(raster_dff(all(~isnan(raster_dff),2),:),1));
        
        XX = [x_all, fliplr(x_all)];
        YY = [dff_mean-dff_ste, fliplr(dff_mean+dff_ste)];

        subplot(3,3,7)
        plot(x_all, dff_mean, 'color', color_465, 'LineWidth', 3); hold on;
        line([0 0], [min(YY), max(YY)], 'Color', [.7 .7 .7], 'LineWidth', 2)

        h = fill(XX, YY, color_465);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('dF/F %', 'FontSize', 12)
        title('465 nm dF/F % (405-subtracted)')
        % c2 = colorbar;

        %% df/f Z-scored
        zscore_dff = zeros(size(raster_dff));
        for dummy_idx = 1:size(raster_dff,1)
            ind = x_all(1,:) < baseline_window(2) & x_all(1,:) > baseline_window(1);
            zb = mean(raster_dff(dummy_idx,ind), 'omitnan'); % baseline period mean
            zsd = std(raster_dff(dummy_idx,ind), 'omitnan'); % baseline period stdev
            zscore_dff(dummy_idx,:)=(raster_dff(dummy_idx,:) - zb)/zsd; % Z score per bin
        end

        % Standard error of the z-score
        zmean = mean(zscore_dff, 'omitnan');
        zerror = std(zscore_dff, 'omitnan')/sqrt(size(zscore_dff(all(~isnan(zscore_dff),2),:),1));

        % Plot heat map
        subplot(3,3,5)

        imagesc(x_all, 1, zscore_dff);
        colormap('jet'); % c1 = colorbar; 
        title('Z-Score Heat Map (405-subtracted)');
        ylabel('Trials', 'FontSize', 12);

        % Fill band values for second subplot. Doing here to scale onset bar
        % correctly
        XX = [x_all, fliplr(x_all)];
        YY = [zmean-zerror, fliplr(zmean+zerror)];

        subplot(3,3,8)
        plot(x_all, zmean, 'color', color_465, 'LineWidth', 3); hold on;
        line([0 0], [min(YY), max(YY)], 'Color', [.7 .7 .7], 'LineWidth', 2)

        h = fill(XX, YY, color_465);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('Z-score', 'FontSize', 12)
        title('465 nm Z-score (405-subtracted)')
        % c2 = colorbar;
        
        %% 465 Z-scored
        zscore_465 = zeros(size(dcOff_465));
        for dummy_idx = 1:size(dcOff_465,1)
            ind = x_all(1,:) < baseline_window(2) & x_all(1,:) > baseline_window(1);
            zb = mean(dcOff_465(dummy_idx,ind)); % baseline period mean
            zsd = std(dcOff_465(dummy_idx,ind)); % baseline period stdev
            zscore_465(dummy_idx,:)=(dcOff_465(dummy_idx,:) - zb)/zsd; % Z score per bin
        end
        % Standard error of the z-score
        zmean = mean(zscore_465, 'omitnan');
        zerror = std(zscore_465, 'omitnan')/sqrt(size(zscore_465(all(~isnan(zscore_465),2),:),1));

        % Plot heat map
        subplot(3,3,6)

        imagesc(x_all, 1, zscore_465);
        colormap('jet'); % c1 = colorbar; 
        title('Z-Score Heat Map (465 only)');
        ylabel('Trials', 'FontSize', 12);

        % Fill band values for second subplot. Doing here to scale onset bar
        % correctly
        XX = [x_all, fliplr(x_all)];
        YY = [zmean-zerror, fliplr(zmean+zerror)];

        subplot(3,3,9)
        plot(x_all, zmean, 'color', color_465, 'LineWidth', 3); hold on;
        line([0 0], [min(YY), max(YY)], 'Color', [.7 .7 .7], 'LineWidth', 2)

        h = fill(XX, YY, color_465);
        set(h, 'facealpha',.25,'edgecolor','none')

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        ylabel('Z-score', 'FontSize', 12)
        title('465 nm Z-score')
        % c2 = colorbar;
        
        output_filename = fullfile(cur_savedir, ['FP_trialPlot_' cur_epoc_name '_code' int2str(cur_epoc_code) '_' cur_epoc_onsetoffset]);
        savefig(f1, [output_filename '.fig'])

        set(gcf, 'PaperPositionMode', 'auto', 'renderer', 'Painters');
        print(gcf, '-painters', '-dpdf', '-r300',  output_filename, '-fillpage')
        
        %% Save dcOffset signals, 465 and 405 fit, df/f, df/f z-score and 465 z-score as individual tables
        datafilepath = split(cur_savedir, filesep);
        subj_id = split(datafilepath{end-1}, '-');
        subj_id = join(subj_id(1:3), "-");
        datafilename = fullfile(cur_savedir, [subj_id{1} '_' ...
            datafilepath{end} '_' cur_epoc_name '_code' ...
            int2str(cur_epoc_code) '_' cur_epoc_onsetoffset]);
        
        % dcOffset signals
        TT = array2table(dcOff_465);
        writetable(TT, [datafilename '_trialData_465sig_dcOffset.csv']);        
        TT = array2table(dcOff_405);
        writetable(TT, [datafilename '_trialData_405sig_dcOffset.csv']);        

        % 465 and 405 fit
        TT = array2table(raster_465);
        writetable(TT, [datafilename '_trialData_465sig_mV.csv']);        
        TT = array2table(raster_405fit);
        writetable(TT, [datafilename '_trialData_405fit_mV.csv']);        

        % dF/F
        TT = array2table(raster_dff);
        writetable(TT, [datafilename '_trialData_465dff.csv']);        
 
        % dF/F z-score
        TT = array2table(zscore_dff);
        writetable(TT, [datafilename '_trialData_465dff_zscore.csv']);        

        % 465 z-score
        TT = array2table(zscore_465);
        writetable(TT, [datafilename '_trialData_465sig_zscore.csv']);    
        
        %% AUC calculation
        % Between -0.5 and 2
        AUC_zscore_dff = trapz(zscore_dff(:, x_all > -0.5 & x_all < 2), 2);
        AUC_zscore_465 = trapz(zscore_465(:, x_all > -0.5 & x_all < 2), 2);
        
        %Compile and save tables

        trial_counter = 1:size(AUC_zscore_dff, 1);
        TT = array2table([trial_counter' AUC_zscore_dff AUC_zscore_465],...
            'VariableNames',{'Trial' 'AUC_zscore_dff' 'AUC_zscore_465'});
        writetable(TT, [datafilename '_AUC.csv']);
    end
    
    tEnd = toc(t0);
    fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
end