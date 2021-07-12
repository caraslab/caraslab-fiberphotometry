function caraslab_getAMDepthData_FPdata(Savedir, sel)

%% Timestamps of interest
% Change this if needed; It will skip non-present ones
% Organized in {NAME, CODE, ONSET/OFFSET}
epocs_names = {{'TTyp', 0, 'onset'}};

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
    fprintf('Getting AM trial data, %s.......\n', cur_path.name)
    
    t0 = tic;
    
    % Load data CSV files
    cur_datafiledir = dir(fullfile(cur_savedir, '*_TTyp_code0_onset_trialData_465dff_zscore.csv'));
    fullpath = fullfile(cur_datafiledir.folder, cur_datafiledir.name);
    dffzscore_df = readtable(fullpath);
    
    cur_datafiledir = dir(fullfile(cur_savedir, '*_TTyp_code0_onset_trialData_465dff.csv'));
    fullpath = fullfile(cur_datafiledir.folder, cur_datafiledir.name);
    dff_df = readtable(fullpath);
    
    cur_datafiledir = dir(fullfile(cur_savedir, ['*' cur_path.name '_TTyp_code0_onset_AUC.csv']));
    fullpath = fullfile(cur_datafiledir.folder, cur_datafiledir.name);
    auc_df = readtable(fullpath);  
    
    % Load behavioral file
    cur_datafiledir = dir(fullfile(cur_savedir, 'CSV files', ['*' cur_path.name '_trialInfo.csv']));
    fullpath = fullfile(cur_datafiledir.folder, cur_datafiledir.name);
    trialInfo_df = readtable(fullpath);      
    trialInfo_AMtrials = trialInfo_df(trialInfo_df.TrialType == 0,:);
    
    % Separate trials by AM depth
    all_AMdepth = sort(unique(trialInfo_AMtrials.AMdepth(trialInfo_AMtrials.Reminder == 0)));
    log_AMdepth = round(make_stim_log(all_AMdepth));
    
    % For each AM depth
    color_resp = [0, 128, 0]/255;
    avg_AUC = [];  % each column is one AMtrial
    ste_AUC = [];  % each column is one AMtrial
    dprime_dffAUC = [];
    dprime_dffzscoreAUC = [];
    
    subplot_counter = 1;
    x_all = linspace(plot_window(1), plot_window(2), size(dffzscore_df, 2));
    close all
    f1 = figure;

    all_resp_ax = [];
    all_heatmap_ax = [];
    for AMdepth_ind=1:length(all_AMdepth)
        AMdepth = all_AMdepth(AMdepth_ind);
        all_resp_ax = [all_resp_ax subplot(4, length(all_AMdepth), subplot_counter)];
        
        % Filter current trials
        cur_ind = trialInfo_AMtrials.AMdepth == AMdepth;
        cur_responses = table2array(dffzscore_df(cur_ind,:));
        
        % Plot mean responses
        mean_resp = mean(cur_responses, 1, 'omitnan');
        ste_resp = std(cur_responses, 1, 'omitnan')/sqrt(size(cur_responses(all(~isnan(cur_responses),2),:),1));
        

        % Create the standard error bands for the 405 signal
        XX = [x_all, fliplr(x_all)];
        YY = [mean_resp + ste_resp, fliplr(mean_resp - ste_resp)];

        % Plot filled standard error bands.
        h = fill(XX, YY, color_resp); hold on;
        set(h, 'facealpha',.25,'edgecolor','none')
        
        plot(x_all, mean_resp, 'color', color_resp, 'LineWidth', 1);

        % Finish up the plot
        axis tight
        xlabel('Time, s','FontSize',12)
        
        if AMdepth_ind ==1
            ylabel('dF/F z-score', 'FontSize', 12)
        end
        
        title(sprintf('AM = %d dB re:100%', round(make_stim_log(AMdepth))))
        
        % Plot heatmaps
        all_heatmap_ax = [all_heatmap_ax subplot(4, length(all_AMdepth), subplot_counter+length(all_AMdepth))];
        imagesc(x_all, 1, cur_responses, [-10 10]);
        colormap('jet'); % c1 = colorbar; 
        title('dF/F Heat Map (405-subtracted)');
        ylabel('Trials', 'FontSize', 12);
        subplot_counter = subplot_counter + 1;
        
        % Gather AUCs to plot later
        cur_aucs = auc_df.AUC_zscore_dff(cur_ind,:);
        avg_AUC = [avg_AUC mean(cur_aucs, 'omitnan')];
        ste_AUC = [ste_AUC std(cur_aucs, 1, 'omitnan')/sqrt(size(cur_aucs(all(~isnan(cur_aucs),2),:), 1))];
        
        % Gather dff_zscores and dffs to get FP d' values
        % dff_zscore
        baseline_AUCs = trapz(cur_responses(:, x_all >= -1 & x_all < 0), 2);
        response_AUCs = trapz(cur_responses(:, x_all >= 0 & x_all < 2), 2);
        dprime = 2*(mean(response_AUCs, 'omitnan') - mean(baseline_AUCs, 'omitnan')) ...
            / (std(response_AUCs, 'omitnan') + std(baseline_AUCs, 'omitnan'));
        dprime_dffzscoreAUC = [dprime_dffzscoreAUC dprime];
        
        % dff
        cur_responses = table2array(dff_df(cur_ind,:));
        baseline_AUCs = trapz(cur_responses(:, x_all >= -1 & x_all < 0), 2);
        response_AUCs = trapz(cur_responses(:, x_all >= 0 & x_all < 2), 2);
        dprime = 2*(mean(response_AUCs, 'omitnan') - mean(baseline_AUCs, 'omitnan')) ...
            / (std(response_AUCs, 'omitnan') + std(baseline_AUCs, 'omitnan'));
        dprime_dffAUC = [dprime_dffAUC dprime];
                
    end
    linkaxes(all_resp_ax,'y')
    linkaxes(all_heatmap_ax,'y')    
    
    % Plot AUCs
    subplot(4, 1, 3)
    hBar = bar(avg_AUC, 'FaceColor', [.8 .8 .8]); hold on;
    
    % Get the x coordinate of the bars
    x = hBar.XEndPoints;

    % Plot the errorbars
    errorbar(x',avg_AUC, ste_AUC,'k','linestyle','none');
    
    % Finish plot
    set(gca,'xticklabel',{log_AMdepth});
    ylabel('AUC', 'FontSize', 12);
    xlabel('AM depth (dB re:100%)', 'FontSize', 12);
    title('Area-under-curve (-0.5 -> 2 s)')
    
    % Plot dprimes
    subplot(4, 1, 4)

    plot(log_AMdepth, dprime_dffzscoreAUC,'-o', 'color', 'k', 'MarkerFaceColor', 'k'); hold on;
    plot([min(log_AMdepth) max(log_AMdepth)], [1 1], 'r:', 'linewidth', 2)
    
    % Finish plot
    ylabel('d''', 'FontSize', 12);
    ylim([min(dprime_dffzscoreAUC) - 0.5 max(dprime_dffzscoreAUC) + 0.5]);
    xlabel('AM depth (dB re:100%)', 'FontSize', 12);
    title('Neural d''')
    set(gcf, 'Position',[100, 100, 2000, 1000])      
    
    % Save fig
    datafilepath = split(cur_savedir, filesep);
    subj_id = split(datafilepath{end-1}, '-');
    subj_id = join(subj_id(1:3), "-");

    output_filename = fullfile(cur_savedir, [subj_id{1} '_' ...
            datafilepath{end} '_AMdepth_responses' ]);

    savefig(f1, [output_filename '.fig'])

    set(gcf, 'PaperPositionMode', 'auto', 'renderer', 'Painters');
    print(gcf, '-painters', '-dpdf', '-r300',  output_filename)
    
    %Compile and save tables

    TT = array2table([log_AMdepth avg_AUC' dprime_dffAUC' dprime_dffzscoreAUC'],...
        'VariableNames',{'Stimulus' 'average_AUC' 'dff_d_prime' 'dffzscore_d_prime'});
    writetable(TT, fullfile(cur_savedir, [subj_id{1} '_' ...
            datafilepath{end} '_neuralDprime.csv']));
    
    tEnd = toc(t0);
    fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
end