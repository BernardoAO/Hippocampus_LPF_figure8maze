%% Td analysis functions
% Bernardo AO

function [vel_grid, ratio_grid, tran, ratio_vel] = ...
    theta_delta_ratio(pideals, lfps, opt)
    % Inputs:
    % pideals: cell of 1d arrays with the pideals object for each session
    % lfps: cell of 1d arrays with the LFP object for each session
    % % opt | extra parametes ( Refer to ITPC_analysis)

    % Outputs:
    % vel_grid: 2d array with the mean velocity for each grid point
    % ratio_grid: 2d array with the mean theta delta ratio
    %   for each grid point
    % tran: down and up states' transition indices.
    
    n_tw_tran = 10;
    thr_tran = 3;
    
    n_sess = numel(lfps);
    vel_vec = [];
    ratio_vec = [];
    x_vec = [];
    y_vec = [];

    for n = 1:n_sess
        lfp = lfps{n};
        pideal = pideals{n};

        log_ratio_pt = get_ratio_pi(pideal, lfp, opt.n_cycles, opt.win);
        vel_vec = [vel_vec; pideal.v];
        ratio_vec = [ratio_vec; log_ratio_pt];

        x_vec = [x_vec; pideal.x];
        y_vec = [y_vec; pideal.y];
    end
    
    tran_u = get_state_transitions(ratio_vec, n_tw_tran, thr_tran);
    tran_d = get_state_transitions(ratio_vec, n_tw_tran, -thr_tran);
    tran = struct('up', tran_u, 'down', tran_d);
    
    ratio_vel = corr_speed_td(vel_vec,ratio_vec,opt);

    % Get grids and plot
    cmap = get_cmap();
    vel_grid = fig8_heatmap(x_vec, y_vec, vel_vec, [], "velocity", ...
        hot, opt);
    ratio_grid = fig8_heatmap(x_vec, y_vec, ratio_vec, tran, "log ratio", ...
        cmap, opt);
    tran = hist_times(x_vec, y_vec, tran, cmap, opt);
end

function conv_result = conv_LFP(lfp, fs, freq_band, n_cycles)
    % Inputs:
    % lfp: 1d array with the LFP signal
    % fs: sampling time_vec
    % freq_band: desired freq_band [start, end]
    % n_cycles: number of cycles for the morlet wavelet
    % 
    % Outputs:
    % conv_result: convolution of the LFP with a morlet wavelett
    
    % Morlet wavelet
    center_freq = round(mean(freq_band));              
    
    w_len = n_cycles / center_freq; 
    t = -w_len/2 : 1/fs : w_len/2;
    s = n_cycles / (2*pi*center_freq); % std of Gaussian
    morlet_wavelet = exp(2*1i*pi*center_freq*t) .* exp(-t.^2/(2*s^2));
    morlet_wavelet = morlet_wavelet / sqrt(sum(abs(morlet_wavelet).^2));
    
    % Time domain
    %conv_result = conv(lfp, morlet_wavelet, 'same');

    % FFT
    delay = floor(length(morlet_wavelet) / 2);
    conv_result = circshift(fftfilt(morlet_wavelet, lfp), -delay);
end

function log_ratio_pt = get_ratio_pi(pideal, lfp, n_cycles, win)
    % Inputs:
    % pideal: pideal object 
    % lfp: LFP object
    % n_cycles: number of cycles for the morlet wavelet
    % win: window for each side to smooth the result

    % Outputs:
    % log_ratio_pt: array with the value of the theta delta
    %   ratio for each pideal value

    theta_band = [6, 12];
    delta_band = [1, 6];

    lfpd = Data(lfp);
    lfpt = EEG.toSecond(lfp, 'Range');
    pt = pideal.t;
    fs = EEG.Fs(lfp);

    tp = abs(conv_LFP(lfpd, fs, theta_band, n_cycles)).^2;
    dp = abs(conv_LFP(lfpd, fs, delta_band, n_cycles)).^2;

    log_ratio = log(tp ./ dp);
    

    log_ratio_pt = zeros(size(pt));
    t = 1;
    T = length(lfpt);
    for i = 1:length(pt)
        while t < T && abs(lfpt(t + 1) - pt(i)) < abs(lfpt(t) - pt(i))
            t = t + 1;
        end
        if t-win > 0
            start_ind = t-win;
        else
            start_ind = 1;
        end
        if t+win < T
            end_ind = t+win;
        else
            end_ind = T;
        end
        log_ratio_pt(i) = mean(log_ratio(start_ind:end_ind));
    end
end

function ratio_down = get_state_transitions(ratio_vec, tw, thr)
    ratio_down = [];
    small_r = ratio_vec < thr;

    for t = 2:length(ratio_vec) - tw
        if small_r(t-1) == 0 && sum(small_r(t:t+tw)) > tw
            ratio_down = [ratio_down, t];
        end
    end
end

function cmap = get_cmap()
    % Color map for ratio
    cn = 256;
    top = linspace(1,0,cn/2)';
    bot = linspace(0,1,cn/2)';
    cmap = [ones(cn/2,1), bot, bot; top, top, ones(cn/2,1)];
end

function z_grid = fig8_heatmap(x_vec, y_vec, z_vec, tran, name, cmap,opt)
    % Gets and plots a heatmap given x,y, and z values. 
    % Plots a scatter as well if given a points index array

    x_edges = -22.5:1:22.5; 
    y_edges = -35:35; 
    x_bins = length(x_edges)-1;
    y_bins = length(y_edges)-1;

    [~, ~, x_idx] = histcounts(x_vec, x_edges);
    [~, ~, y_idx] = histcounts(y_vec, y_edges);

    valid = x_idx > 0 & y_idx > 0;
    x_idx = x_idx(valid);
    y_idx = y_idx(valid);
    z_vec = z_vec(valid);

    linear_idx = sub2ind([y_bins, x_bins], y_idx, x_idx);
    z_sum = accumarray(linear_idx, z_vec, [y_bins*x_bins, 1], @sum, NaN);
    z_count = accumarray(linear_idx, 1, [y_bins*x_bins, 1], @sum, NaN);
    z_avg = z_sum ./ z_count;

    z_grid = reshape(z_avg, [y_bins, x_bins]);

    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

    %% Plot
    if ~isempty(tran) % up and down scatter
        s_names = fieldnames(tran);
        z_grid(isnan(z_grid)) = 0;
        for s = 1:length(s_names)

            states = s_names{s};
            points = tran.(states);

            figure;
            imagesc(x_centers, y_centers, z_grid);
            axis xy;
            axis off
            lp = length(points);
            x_tran = x_vec(points) + 0.1*randn(lp,1);
            y_tran = y_vec(points) + 0.1*randn(lp,1);
        
            hold on
            h = scatter(x_tran,y_tran,"filled");
            h.MarkerFaceAlpha = 0.2;
            h.MarkerFaceColor = 'black';
            if strcmp(states, 'up')
                h.Marker = "^";
            elseif strcmp(states, 'down')
                h.Marker = "v";
            end
            hold off
            t_name = opt.fig_name + " Average " + name + " heatmap " + ...
                states + " states";
            title(t_name);
            colormap(flipud(cmap))
            clim([-4,4]);
            a = colorbar;
            a.Label.String = name;
            saveas(gcf, fullfile(opt.save_path, t_name + opt.ext));
        end

    else % no scatter
        figure;
        imagesc(x_centers, y_centers, z_grid);
        axis xy;
        axis off
        t_name = opt.fig_name + " Average " + name + " heatmap";
        title(t_name);
        colormap(flipud(cmap))
        a = colorbar;
        a.Label.String = name;

        if name == "theta power" || name == "delta power"
            m = mean(z_grid,"all","omitmissing");
            s = std(z_grid, [],"all","omitmissing");
            clim([0,m+3*s]);
        end
        saveas(gcf, fullfile(opt.save_path, t_name + opt.ext));
    end
end

function tran = hist_times(x_vec, y_vec, tran, cmap,opt)
    x_edges = -22.5:1:22.5; 
    y_edges = -35:35;
    
    s_names = fieldnames(tran);
    for s = 1:length(s_names)
        states = s_names{s};
        points = tran.(states);
        pointsx = x_vec(points);
        pointsy = y_vec(points);
        [t_counts,~,~] = histcounts2(pointsx, pointsy,x_edges, y_edges);
        
        np = 1000;
    
        total_counts = zeros(length(x_edges)-1,length(y_edges)-1,np);
        for p = 1:np
            r = randperm(length(x_vec),length(points));
            pointsx_r = x_vec(r);
            pointsy_r = y_vec(r);
            [counts,~,~] = histcounts2(pointsx_r, pointsy_r,x_edges, y_edges);
            total_counts(:,:,p) = counts;
        end
    
        m = mean(total_counts,3);
        st = std(total_counts,[],3);
        t_norm_counts = (t_counts - m);
        c = (t_norm_counts ~= 0); 
        t_norm_counts(c) = t_norm_counts(c) ./ (st(c) + 1e-3);
        t_norm_counts = t_norm_counts.';
        tran.(states) = {tran.(states), t_norm_counts};

        x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
        y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;
    
        % Plot
        figure;
        imagesc(x_centers, y_centers, t_norm_counts);
        axis xy;
        axis off
        colormap(flipud(cmap));
        clim([-3,3]);
        a = colorbar;
        a.Label.String = 'Count z-score';
        t_name = opt.fig_name + " Normalized histogram " + ...
                    states + " states";
        title(t_name);
        saveas(gcf, fullfile(opt.save_path, t_name + opt.ext));
    end

    %figure;
    %histogram2(pointsx, pointsy,x_edges, y_edges)
    %colormap(flipud(hot))
    %grid off
    %xticks([])
    %ticks([])
end

function r_v = corr_speed_td(v,r,opt)
    % Plots a scatter and gets the correlation of a velocity vector and a
    % ratio vector (v,r)    
    v_edges = opt.v_edges;
    v_bins = (v_edges(1:end-1) + v_edges(2:end)) / 2;
    r_v = zeros(length(v_bins),1);
    
    [~, ~, bin_indx] = histcounts(v, v_edges);
    for i = 1:length(v_bins)
        r_v(i) = mean(r(bin_indx == i),'omitnan');
    end

end

function plot_example_transitions(pideal, log_ratio_pt, ratio_down, thr)
    figure;
    tv = pideal.t - pideal.t(1);
    plot(tv,log_ratio_pt,Color="#ae2012")
    for t = ratio_down
        xline(tv(t),LineStyle="--",Color="#0a9396")
    end
    yline(thr,LineStyle=":",Color="#343a40")
    xlim([0,50])
    xlabel("time [s]")
    ylabel("log theta/delta")
    box off
end

% Old code
%{

 lfp_raw = zeros(n_trials, len_t);
            lfp_raw(ct, ti) = lfp(align_indx + t);
        r = log(theta_power(trial,:) ./ delta_power(trial,:));
        figure();
        yyaxis left
        plot(time_window_vec, r,color="k")
        yyaxis right
        plot(time_window_vec, theta_power(trial,:),color="blue")
        hold on
        plot(time_window_vec, delta_power(trial,:),color="red")
        hold off

first_t = {};
    s = 1;
    for trial = 1:n_trials
        basal_r = min(ratio(trial,time_window_vec < 0));
        t_small_td = time_window_vec(ratio(trial,:) < basal_r);
        if ~isempty(t_small_td)
            first_t{s} = t_small_td(1);
            s = s + 1;
        end
    end

function [ratio_mean, min_time, max_time] = theta_delta_ratio(lfps, ...
    align_indices, time_window, opt)
    % Inputs:
    % Same as ITPC_analysis

    save_path = "W:\Lorena\Analysis_scripts\Bernardo_code\plots";
    ext = ".png";
    plot_ts = "max";
    theta_band = [6, 12];
    delta_band = [1, 6];

    time_window_ind = int32(time_window(1)*opt.fs: opt.fs*opt.step : time_window(2)*opt.fs);
    time_window_vec = time_window(1):opt.step:time_window(2);

    len_t = length(time_window_ind);
    n_sessions = numel(lfps);
    n_trials = calc_n_trials(align_indices);

    theta_power = zeros(n_trials, len_t); 
    delta_power = zeros(n_trials, len_t); 
    

    %% delta theta power
    trial = 1;
    for n = 1:n_sessions
        lfp = lfps{n};
        align_indx = int32(align_indices{n});
        n_t = length(align_indx);
        ct = trial:trial + n_t - 1;

        
        tp = abs(conv_LFP(lfp, opt.fs, theta_band, 10)).^2;
        dp = abs(conv_LFP(lfp, opt.fs, delta_band, 10)).^2;

        ti = 0;
        for t = time_window_ind
            ti = ti + 1;
            theta_power(ct, ti) = tp(align_indx + t);
            delta_power(ct, ti) = dp(align_indx + t);
        end

        trial = trial + n_t;
    end
    t_mean = mean(theta_power, 1,"omitmissing");
    d_mean = mean(delta_power, 1,"omitmissing");
    ratio_mean = log(t_mean ./ d_mean);

    % find the time of the min ratio
    ratio = log(theta_power ./ delta_power);
    [~,min_ind] = min(ratio,[],2);
    min_time = time_window_vec(min_ind);
    
    [~,max_ind] = max(delta_power,[],2);
    max_time = time_window_vec(max_ind);

    %% Plot
    figure('Name','Theta delta ratio' + opt.band_name);
    subplot(2,1,1)
    p1 = plot(time_window_vec, t_mean, color="blue");
    hold on
    p2 = plot(time_window_vec, d_mean, color="red");
    xline(0,LineStyle=":",Color="#343a40")
    legend([p1 p2], {"theta","delta"})
    xlim(time_window)
    ylim([0,0.5e-4])
    hold off    
    box off
    ylabel("Power")

    subplot(2,1,2)
    plot(time_window_vec, ratio_mean, color=opt.color)
    hold on
    xline(0,LineStyle=":",Color="#343a40")
    if plot_ts == "max"
        ts = max_time;
    else
        ts = min_time;
    end
    for t = ts 
        xline(t,LineStyle=":",Color="green")
    end
    hold off
    xlabel("time [s]")
    xlim(time_window)
    %ylim([-2 2])
    ylabel("log ratio")
    box off
    t_name = opt.band_name + ' delta ratio ' + opt.align_name;
    
    sgtitle(t_name)
    saveas(gcf, fullfile(save_path, t_name + ext));
end

% plot session
figure;
yyaxis left
plot(log_ratio_pt, Color="blue")
hold on
yyaxis right
plot(pideal.v,color="red")
hold off
legend("log_ratio","v")

figure;
h = scatter(pideal.v,log_ratio_pt, Color=opt.color);
h.MarkerFaceAlpha = 0.4;
ylim([0,1e-4])

error("done")

%}
