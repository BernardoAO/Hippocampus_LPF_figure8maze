%%%% TD speed Analysis
%%% Bernardo AO, adapted from ITPC_analysis.m
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
delay_t = {'d10'}; % sessInfo(i).sessDirs

freq_band = [6, 12]; % [6, 12] Theta, [12, 20] Beta
fig_name = delay_t{1} + " Theta";

show_fig = true;

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
    pideals = {};
    
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

            % Speed


            % Collect data
            lfps{s} = lfp;
            pideals{s} = pideal;

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
    if td_ratio
    [ratio_mean, min_time, max_time] = theta_delta_ratio(lfpds, ...
        region_ends, time_window, opt_itpc);
    results(a).ratio_mean = ratio_mean;
    results(a).min_time = min_time;
    results(a).max_time = max_time;
    else
    try
        [band_itpc, band_itpc_t, n_t] = ITPC_analysis(lfpds, region_ends, ...
            time_window, opt_itpc);

        results(a).n_t = n_t;
        if n_t > 20 % metric only stable for n > 40
            results(a).band_itpc = band_itpc;
            results(a).band_itpc_t = band_itpc_t;
        end
    catch
        problematicSessions = [problematicSessions; animal];
    end
    end
    %% Save
    results(a).animal = animal;
    results(a).genotype = sessInfo.genotype;
    
    r = a/length(animal_numbers);
    waitbar(r,wb,"Running " + num2str(r*100,2) + "%");
    a = a + 1;
end
delete(wb);

%% Results

if td_ratio
    %save('td_results.mat','results')
    clc
else
    genotypesc = categorical(genotypes);
    genotypesc = reordercats(genotypesc, genotypes);
    plot_itpc(results, genotypesc, time_window, colors, ...
        opt_itpc.align_name + " " + delay_t{1})
end