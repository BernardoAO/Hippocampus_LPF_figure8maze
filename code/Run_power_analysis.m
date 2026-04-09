%%%% TD speed Analysis
%%% Bernardo AO, adapted from ITPC_analysis.m
% Data needed: sessionInfo.mat, CSC(n).ncs
% Script dependencies: ITPC_analysis.m,  plotTrialPerCell, trialInfo2simple, readCRTsd
% 
addpath(genpath('W:\LABanalysis\SilviaProjectCode\AnalysisPerTrial\RunAnalysis'))
addpath("W:\Bernardo\Hippocampus_LPF_figure8maze")

disp('Starting ITPC Analysis')
clc
% 'W:\Lorena\Analysis_scripts\DataOrganization'
%% Define which animal / days to run 

All_sessInfo = load(fullfile('W:\LABanalysis\SilviaProjectCode\+Figure8DataOrganization', ...
    'sessionInfo.mat')).sessInfo;

animal_numbers = unique([All_sessInfo.animal]);
all_rois = {'return';'delay';'stem';'choice';'reward'};
genotypes = ["Control", "CA1-APP"];
colors = [0.14,0.85,0.71; 0.85,0.14,0.28];
sp = "W:\Bernardo\plots";

% Parameters 
win = 640; % window to smooth the ratio around win / fs = +- 20 ms
v_edges = [0, logspace(-0.6,1.2,12)]; % velocity edges for binning
band = [5, 12]; % theta
smooth_win = 2; % seconds
delay_t = {'d0'}; %sessInfo(i).sessDirs

% Outputs
problematicSessions = {};
results = struct([]);

wb = waitbar(0,'Running animal loop 0%');
a = 1;
animal_numbers = [58];%[9,12,14,19,23,24,27,28,30,38,39,113,114,115,116];
for animal = animal_numbers
    sessions = {All_sessInfo.animal};
    s_indices = find(cellfun(@(sessions) ~isempty(sessions) && ...
        sessions==animal, sessions));
    
    lfps = {};
    praws = {};
    pideals = {};
    
    s = 0;
    for iii = s_indices

        sessInfo = All_sessInfo(iii);
        sessInfo.mainDir = strrep(sessInfo.mainDir, 'W:\Silvia', 'W:\Silvia\UCSD');
        ptempl = parsingtemplate('fig8mouse:rdscr');
                
        i = 1;
        try
            % Liniarize path
            [trialInfo, parsingInfo, pathData, pathDataIdeal, ...
                pathDataLin] = plotTrialPerCell.loadInfo(sessInfo, i);

            s = s + 1;
            for block_c = delay_t 
                block = block_c{1};
                blockDir = fullfile(sessInfo(i).mainDir, block);
                LFP_dir = fullfile(blockDir,'LFP'); if ...
                    ~exist(LFP_dir,'dir'), mkdir(LFP_dir),end
    
                tinfo = trialInfo.(block);
                praw = pathData.(block);
                pideal = pathDataIdeal.(block);
   
                %% Collect LFPs
                channelInLayer = sprintf('CSC%d.ncs', sessInfo.cellLayerChann); % Picks the channel that is in the layer
                lfp = readCRTsd(fullfile(blockDir, channelInLayer));
    
                lfps{s} = lfp;
                praws{s} = praws;
                pideals{s} = pideal;

            end
        catch
            problematicSessions = [problematicSessions; sessInfo(i).mainDir];
            continue;
        end
        clc
    end
    
    %% Run power analysis
    fprintf('\nRunning power analysis\n');
    opt = struct(...
                'fs', EEG.Fs(lfp), ...
                'band', band, ...
                'color', colors(sessInfo.genotype,:), ...
                'fig_name', int2str(animal), ...
                'smooth_win', smooth_win, ...
                'v_edges', v_edges, ...
                'save_path', sp, ...
                'ext', ".png");
    
    
    %try
    [vel_grid, power_grid, ratio_vel] = ...
        power_analysis(pideals, lfps, opt);

    results(a).vel_grid = vel_grid;
    results(a).power_grid = power_grid;
    results(a).ratio_vel = ratio_vel;

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
clc
delete(wb);
error "end"


%%
plot_mean_grid(results,"power_grid", "theta power", opt_td)
for g = 1:2
    r = results([results.genotype] == g);
    plot_mean_grid(r,"power_grid", "theta power " + genotypes(g), opt_td)
end
%% Results

plot_mean_grid(results,"vel_grid", "velocity", opt_td)
plot_mean_grid(results,"ratio_grid", "norm log ratio ", opt_td)
plot_mean_grid(results,"norm_tran_down", "norm transition down ", opt_td)
plot_mean_grid(results,"norm_tran_up", "norm transition up ", opt_td)

%% Groups

for g = 1:2
    r = results([results.genotype] == g);
    plot_mean_grid(r,"vel_grid", "velocity " + genotypes(g), opt_td)
    plot_mean_grid(r,"ratio_grid", "norm log ratio " + genotypes(g), opt_td)
    plot_mean_grid(r,"norm_tran_down", "norm transition down " ...
       + genotypes(g), opt_td)
    plot_mean_grid(r,"norm_tran_up", "norm transition up " ...
       + genotypes(g), opt_td)
end

%% Groups

plot_ratio_velocity(results, colors, genotypes, opt_td)

function plot_mean_grid(results, gname, name, opt)
    
    x_edges = -23:23;
    y_edges = -35:35;
    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
    
    g_shape = size(results(1).(gname));
    n_animals = numel(results);
    all_grid = zeros(g_shape(1), g_shape(2), n_animals);

    for a = 1:n_animals
        if contains(name,"theta power") || contains(name,"delta power")
            grid = results(a).(gname);
            grid(grid==0) = nan;
            m = mean(grid,"all","omitmissing");
            s = std(grid, [],"all","omitmissing");
            all_grid(:,:,a) = (grid - m) / s;
        else
            all_grid(:,:,a) = results(a).(gname);
        end
    end
    mean_grid = mean(all_grid,3);

    % Plot
    figure;
    imagesc(x_centers, y_centers, mean_grid);
    axis xy; 
    t_name = " Average " + name + " heatmap";
    title(t_name);
    axis off
    if gname == "vel_grid" || contains(name,"theta power") || contains(name,"delta power")
        colormap(flipud(hot))
        if gname == "vel_grid"
            clim([0,14]);
        else
            clim([-1,3]);
        end
    else
        cmap = get_cmap();
        colormap(flipud(cmap))
        clim([-2,2]);
    end
    a = colorbar;
    a.Label.String = name;
    saveas(gcf, fullfile(opt.save_path, t_name + "all.pdf"));
end

function plot_ratio_velocity(results, colors, genotypes, opt)
    v_edges = opt.v_edges;
    v_bins = (v_edges(1:end-1) + v_edges(2:end)) / 2;
    
    figure;
    for g = 1:2
        r = results([results.genotype] == g);
        m = mean([r.ratio_vel], 2);
        err = std([r.ratio_vel], [], 2) / sqrt(size([r.ratio_vel],2));
        hold on
        e = errorbar(v_bins, m, err);
        e.Color = colors(g,:);
        e.Marker = "o";
    end
    hold off
    legend(genotypes)
    xlabel("velocity")
    ylabel("log td ratio")
    t_name = " Average td ratio";
    title(t_name);
    box off

    saveas(gcf, fullfile(opt.save_path, t_name + "all.pdf"));
end

function cmap = get_cmap()
    % Color map for ratio
    cn = 256;
    top = linspace(1,0,cn/2)';
    bot = linspace(0,1,cn/2)';
    cmap = [ones(cn/2,1), bot, bot; top, top, ones(cn/2,1)];
end