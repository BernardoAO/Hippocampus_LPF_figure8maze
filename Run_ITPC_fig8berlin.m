%%%% ITPC Analysis
%%% Bernardo AO, adapted from ...
% Data needed: sessionInfo.mat, CSC(n).ncs
% Script dependencies: ITPC_analysis.m,  plotTrialPerCell, trialInfo2simple, readCRTsd
% 
addpath(genpath('W:\LABanalysis\SilviaProjectCode\AnalysisPerTrial\RunAnalysis'))
addpath("W:\Lorena\Analysis_scripts\Bernardo_code\")

disp('Starting ITPC Analysis')
clc

%% Define which animal / days to run 

load(fullfile('W:\Lorena\Analysis_scripts\DataOrganization', 'sessionInfo.mat'));
All_sessInfo = sessInfo; clear sessInfo
animal_numbers = unique([All_sessInfo.animal]);
all_rois = {'return';'delay';'stem';'choice';'reward'};
genotypes = ["Control", "CA1-APP"];
colors = [0.14,0.85,0.71; 0.85,0.14,0.28];

% Parameters 
region = "stem"; % align to the region end
time_window = [-1.5, 1.5]; % for ITPC around the onset

delay_t = {'d10'}; % sessInfo(i).sessDirs

freq_band = [6, 12]; % [6, 12] Theta, [12, 20] Beta
fig_name = delay_t{1} + " Theta";

show_fig = false;

% Outputs
problematicSessions = {};
results = struct([]);

wb = waitbar(0,'Running animal loop 0%');
a = 1;
%animal = 1433; sessInfo = All_sessInfo(1); block_c = sessInfo(i).sessDirs 
for animal = animal_numbers
    sessions = {All_sessInfo.animal};
    s_indices = find(cellfun(@(sessions) ~isempty(sessions) && ...
        sessions==animal, sessions));
    
    lfpds = {}; 
    region_ends = {};
    
    s = 0;
    for iii = s_indices

        sessInfo = All_sessInfo(iii);        
        ptempl = parsingtemplate('fig8mouse:rdscr');
                
        i = 1;
        try
            % Liniarize path
            [trialInfo, parsingInfo, pathData, pathDataIdeal, ...
                pathDataLin] = plotTrialPerCell.loadInfo(sessInfo, i);

            s = s + 1;
            for block_c = delay_t % sessInfo(i).sessDirs
            block = block_c{1};
            blockDir = fullfile(sessInfo(i).mainDir, block);
            LFP_dir = fullfile(blockDir,'LFP'); if ...
                ~exist(LFP_dir,'dir'), mkdir(LFP_dir),end

            tinfo = trialInfo.(block);
            pideal = pathDataIdeal.(block);

            tridx = trialInfo2simple.trial_startend_ind(tinfo);
            tridx = tridx(~tinfo.degen, :);
%pinfo = parsingInfo.(block);%plin = pathDataLin.(block);%pdata = pathData.(block);%trtype = plotTrialPerCell.categtable(tinfo);
            
            %% Processing LFPs
            channelInLayer = sprintf('CSC%d.ncs', sessInfo.cellLayerChann); % Picks the channel that is in the layer
            lfp = readCRTsd(fullfile(blockDir, channelInLayer));
            roixy = ptempl(strcmp({ptempl.zone}, region));

            [lfpx_region, ~, ~, ~] = plotTrialPerCell.extractEEGEpoch(pideal, ...
                tridx, roixy, lfp, [0, 0]); %lfpidx_region, lfptime_region, tdomain_bin

            % Collect data for ITPC 
            lfpd = Data(lfp);
            lfpd = lfpd - mean(lfpd);
            lfpds{s} = lfpd;

            region_end = arrayfun(@(i) lfpx_region(i).idx(2), ...
                1:length(lfpx_region));
            region_ends{s} = region_end;

            end
        catch
            problematicSessions = [problematicSessions; sessInfo(i).mainDir];
            %errordlg(sprintf('Error in linearizing path of %s', num2str(iii)));
            continue;
        end
        clc
    end
    
    %% Run ITPC
    fprintf('\nRunning ITPC\n');
    opt_itpc = struct(...
                'fs', lfpx_region(1).Fs, ...
                'step', 0.01, ...
                'freq_band', freq_band, ...
                'n_cycles', 6, ...
                'band_name', int2str(animal) + fig_name, ...
                'align_name', region + " end", ...
                'color', colors(sessInfo.genotype,:), ...
                'show_fig', show_fig);

    %try
        [band_itpc, band_itpc_t, n_t] = ITPC_analysis(lfpds, region_ends, ...
            time_window, opt_itpc);

        results(a).n_t = n_t;
        if n_t > 20 % metric only stable for n > 40
            results(a).band_itpc = band_itpc;
            results(a).band_itpc_t = band_itpc_t;
        end
    %catch
    %    problematicSessions = [problematicSessions; animal];
    %end

    %% Save
    results(a).animal = animal;
    results(a).genotype = sessInfo.genotype;
    
    r = a/length(animal_numbers);
    waitbar(r,wb,"Running " + num2str(r*100,2) + "%");
    a = a + 1;
end
delete(wb);

%% Results

genotypesc = categorical(genotypes);
genotypesc = reordercats(genotypesc, genotypes);
plot_itpc(results, genotypesc, time_window, colors, ...
    opt_itpc.align_name + " " + delay_t{1})


function plot_itpc(results, genotypes, time_window, colors, name)
    %% get values for plotting
    save_path = "W:\Lorena\Analysis_scripts\Bernardo_code\plots";
    gv = [results.genotype];
    tv = linspace(time_window(1), time_window(2), ...
        length(results(1).band_itpc_t));
    itpc_genotype = zeros(length(genotypes),1);
    mean_t = zeros(length(genotypes),length(tv));
    std_t = zeros(length(genotypes),length(tv));
    y = {};
    y_t = {};

    figure('Name','All bar');
    for g = 1:length(genotypes)
        y{g} = [results(gv == g).band_itpc];
        y_t{g} = [results(gv == g).band_itpc_t];
        itpc_genotype(g) = mean(y{g},"omitmissing");
        mean_t(g,:) = mean(y_t{g},2,"omitmissing");
        std_t(g,:) = std(y_t{g},0,2,"omitmissing") / sqrt(size(y_t{g},2));
    end
    
    %% bar plot
    b = bar(genotypes,itpc_genotype, FaceColor="flat");
    b.CData(:,:) = colors;
    hold on
    for g = 1:length(genotypes)
        scatter(genotypes(g), y{g}, "filled", MarkerFaceColor='k')
    end
    
    hold off
    ylabel("ITPC")
    t_name = "Theta " + " ITPC at " + name;
    title(t_name)
    saveas(gcf, fullfile(save_path, t_name + '.pdf'));

    %% time plot
    figure('Name','All time');
    h = [];
    xline(0,LineStyle=":",Color="#343a40")
    hold on;
    yline(0,LineStyle="--")
    for g = 1:length(genotypes)
        upper = mean_t(g,:) + std_t(g,:);
        lower = mean_t(g,:) - std_t(g,:);
        
        fill([tv, fliplr(tv)], [upper, fliplr(lower)], ...
             colors(g,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');     
        
        p = plot(tv, mean_t(g,:), 'Color', colors(g,:), 'LineWidth', 2);
        h(end+1) = p;        
    end
    legend(h, genotypes)
    ylabel("ITPC")
    ylim([-2,3])
    xlabel("time [t]")
    t_name = "Theta " + " ITPC around " + name;
    title(t_name)
    saveas(gcf, fullfile(save_path, t_name + '.pdf'));
end

%Old code
%{
            for trial = 1:length(tridx)
                %ROI_collect_duration = [];all_collect_duration= [];
                for region = 1:5
                    roixy = ptempl(strcmp({ptempl.zone}, all_rois{region}));
                    %start_trial = tridx(trial,1); stop_trial = tridx(trial,1);
                    % trial_member = ismember(bin_member,[tridx(trial,1):tridx(trial,2)])
                    trialx = [tridx(trial,1) tridx(trial,2)];
                    [eegx_region, eegidx_region, eegtime_region, tdomain_bin] = plotTrialPerCell.extractEEGEpoch(pideal, trialx, roixy, eeg, [0, 0]);  %isolate eeg from this area
                    
                    
                    error("done")
                    %remove trials if an area was completely skipped/
                    %tracking missing
                    if ~isstruct(eegx_region)  % animal not tracked in this part of the maze
                        errordlg(sprintf('Could extract EEG in region %d, lap %d animal/i %d, will erase this lap', region, trial, iii));
                        ROI_duration_time(trial,region) = NaN;
                        continue
                    else % regular case
                        times = eegtime_region(2) - eegtime_region(1);
                        all_regions_eegx{trial,region} = eegx_region;
                        all_regions_eegxidx{trial,region} = eegidx_region;
                        all_regions_times{trial,region} = times;
                        all_regions_tdomainbin{trial,region} = tdomain_bin;
                        %Feb 19 2018 added velocity
                        velocity{trial,region} = pdata.v(tdomain_bin.startidx:tdomain_bin.endidx);
                        vid_ts{trial,region} = pdata.t(tdomain_bin.startidx:tdomain_bin.endidx);
                        
                        %% Used for spectogramm
                        fprintf('.');
                        [SPG, t, f, bandSpecgramFun] = specgramwwd(eegx_region.d,eegx_region.Fs, 2, 300,opt.wavelet);
                        spect{trial,region} = SPG;
                        frequ{trial,region} = f;
                        ROI_duration_time(trial,region) = eegtime_region(2) - eegtime_region(1);
                        
                        
                    end
                    %filename = fullfile(savefolder,sprintf('EEG_%s_Animal_i%03d_%s_Region%d.mat', group, iii,block{1},region));
                    %save(filename, 'eegx_region');
                    %
                end % For this ROI
            end
            fprintf('.\n\n');
            %% Now delete trials that have missing data
            delete_row = isnan(mean(ROI_duration_time,2));
            all_regions_eegx(find(delete_row),:) = [];
            all_regions_eegxidx(find(delete_row),:) = [];
            all_regions_times(find(delete_row),:) = [];
            all_regions_tdomainbin(find(delete_row),:) = [];
            velocity(find(delete_row),:) = [];
            vid_ts(find(delete_row),:) = [];
            spect(find(delete_row),:) = [];
            frequ(find(delete_row),:) = [];
            
            successful = tinfo.success(~delete_row);%(1:size(ROI_duration_time,1));
            ROI_duration_time(find(delete_row),:) = [];
            
            %% speed_amp_speed_frequ
            % ->       j = batch('speed_amp_speed_frequ_Regr','Pool',4,'CaptureDiary',true);
            % ->      wait(j);   % Wait for the job to finish
            %      diary(j)   % Display the diary
            %     load(j)    % Load job workspace data into client workspace
            
            %% Save
            filename = fullfile(LFP_dir,sprintf('%s_Animal_i%03d_%s.mat', group, iii,block{1}));
            save(filename,'sessInfo','ROI_duration_time','successful','spect','frequ','all_regions_eegxidx','all_regions_times','all_regions_tdomainbin','velocity', 'vid_ts');
            %% Save info in metafile for analysis
            metafilename=fullfile('W:\Lorena\Analysis_scripts\DataOrganization\metafile_RegionalEEGs.mat');
            if exist(metafilename,'file'), load(metafilename);end
            datainfo.group = group;
            datainfo.animal = sessInfo.animal;datainfo.age = sessInfo.age; datainfo.session = sessInfo.session;     datainfo.filename = filename;
            metainfo.(block{1})(iii) = datainfo;
            %metainfo.(block{1})(iii).z_score_red_file = filename_reduced_z_score;
            %metainfo.(block{1})(iii).z_score_red_file = filename_z_score;
            
            save(metafilename,'metainfo');
            
            
        end % for block
        fprintf('Processing time for Animal No.=%d: ',iii);
        toc(timer_anim);
        %catch
        %failed_iii = [failed_iii,iii];
        %errordlg(sprintf('Could not process i %d Animal %d Day %d', iii,sessInfo.animal,sessInfo.day));
        %continue
        %end %try
    end 
    
end  % end of cycling through days

elapsed_min = round(etime(clock,timer_prog)/60);
msgbox(sprintf('Total runtime: %d min\n',elapsed_min),'Runtime');
%}