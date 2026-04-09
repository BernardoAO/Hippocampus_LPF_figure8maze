%% Td analysis functions
% Bernardo AO

function [vel_grid, power_grid, ratio_vel] = ...
    power_analysis(pideals, lfps, opt)
    % Inputs:
    % pideals: cell of 1d arrays with the pideals object for each session
    % lfps: cell of 1d arrays with the LFP object for each session
    % % opt | extra parametes (Refer to Run_)

    % Outputs:
    % vel_grid: 2d array with the mean velocity for each grid point
    % power_grid: 2d array with the mean power for the given band width
    %             for each grid point
    
    n_sess = numel(lfps);
    vel_vec = [];
    power_vec = [];
    x_vec = [];
    y_vec = [];
    d = designfilt('bandpassiir', ...
    'FilterOrder',4, ...
    'HalfPowerFrequency1',opt.band(1), ...
    'HalfPowerFrequency2',opt.band(2), ...
    'SampleRate',opt.fs);

    for n = 1:n_sess
        lfp = lfps{n};
        lfpd = Data(lfp);
        pideal = pideals{n};
        
        vel_vec = [vel_vec; pideal.v];

        [power_pt, ~] = LFP_Hilbert(lfpd, d, opt);
        
        power_vec = [power_vec; power_pt];

        x_vec = [x_vec; pideal.x];
        y_vec = [y_vec; pideal.y];
    end
    
    ratio_vel = corr_speed_td(vel_vec, power_vec, opt);

    % Get grids and plot
    cmap = get_cmap();

    vel_grid = fig8_heatmap(x_vec, y_vec, vel_vec, [], "velocity", ...
        hot, opt);
    %vel_lin = fig8_linear(vel_grid, "speed");

    power_grid = fig8_heatmap(x_vec, y_vec, power_vec, [], "power", ...
        hot, opt);
    %power_lin = fig8_linear(power_grid, "theta power");
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

function power = power_LFP_multitaper(lfp, fs, freq_band, TW, K)
    % Inputs:
    % lfp: 1d array with the LFP signal
    % fs: sampling time_vec
    % freq_band: desired freq_band [start, end]
    % TW: time-bandwidth produc
    % K: number of tapers
    % 
    % Outputs:
    % power: power of the LFP 

    T = length(lfp);

    % Generate DPSS tapers
    [tapers, ~] = dpss(T, TW, K);

    % FFT parameters
    freqs = (0:T-1)*(fs/T);
    freq_idx = freqs >= freq_band(1) & freqs <= freq_band(2);

    mt_power = zeros(K, T);

    for k = 1:K
        tapered_signal = lfp(:) .* tapers(:,k);
        X = fft(tapered_signal);

        % Band-limited analytic signal
        X_band = zeros(size(X));
        X_band(freq_idx) = X(freq_idx);
        X_band(end-freq_idx+1) = conj(X(freq_idx)); % mirror

        mt_power(k,:) = abs(ifft(X_band)).^2;
    end

    % Average across tapers
    power = mean(mt_power, 1);
end

function [power_sm, phase] = LFP_Hilbert(lfp, d, opt)
    % Inputs:
    % lfp: 1d array with the LFP signal
    % b, a: butter filter 
    % % opt | extra parametes (band, fs, smooth_win)

    % Outputs:
    % power_sm: power smoothed in time of the desired band
    % phase: phase of the desired band

    band_lfp = filtfilt(d, lfp);

    hilbert_lfp = hilbert(band_lfp);
    
    power = abs(hilbert_lfp).^2;
    power_sm = movmean(power, round(opt.smooth_win * opt.fs));

    phase = angle(hilbert_lfp);
end

function power_pt = get_power(pideal, lfp, opt)
    % Inputs:
    % pideal: pideal object 
    % lfp: LFP object
    % % opt | extra parametes (Refer to Run_)

    % Outputs:
    % power_pt: array with the value of the power for each pideal value

    lfpd = Data(lfp);
    lfpt = EEG.toSecond(lfp, 'Range');
    pt = pideal.t;
    fs = EEG.Fs(lfp);

    power = power_LFP_multitaper(lfpd, fs, opt.freq_band, opt.TW, opt.K);
    
    power_pt = zeros(size(pt));
    t = 1;
    T = length(lfpt);

    for i = 1:length(pt)
        while t < T && abs(lfpt(t + 1) - pt(i)) < abs(lfpt(t) - pt(i))
            t = t + 1;
        end
        if t-win > 0
            start_ind = t - opt.win;
        else
            start_ind = 1;
        end
        if t+win < T
            end_ind = t + opt.win;
        else
            end_ind = T;
        end
        power_pt(i) = mean(power(start_ind:end_ind));
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

function linear_grid = fig8_linear(grid, name)
    armx_left = 1;
    armx_center = 23;
    armx_right = 45;
    y_top = 1;
    y_bottom = 70;
    returny_start = 35;
    x_ticks = [0, 59, 85, 127, 155];
    x_labels = ["return", "delay", "stem", "choice", "reward"];
    linear_grid = [];
    
    % Return
    linear_grid = [linear_grid, mean([grid(returny_start:y_bottom-1, armx_left) ...
                                      grid(returny_start:y_bottom-1, armx_right)]')];
    linear_grid = [linear_grid, mean([grid(y_bottom-1, armx_left+1) ...
                                      grid(y_bottom-1, armx_right-1)])];
    linear_grid = [linear_grid, mean([grid(y_bottom, armx_left+1:armx_center); ...
                                      grid(y_bottom, flip(armx_center:armx_right-1))])];
    
    % Center arm
    linear_grid = [linear_grid, grid(flip(y_top+1:y_bottom-1), armx_center)'];
    
    % Top arm
    linear_grid = [linear_grid, mean([grid(y_top, armx_center:armx_right-1); ...
                                      grid(y_top, flip(armx_left+1:armx_center))])];
    
    % Reward
    linear_grid = [linear_grid, mean([grid(y_top+1, armx_left+1) ...
                                      grid(y_top+1, armx_right-1)])];
    linear_grid = [linear_grid, mean([grid(y_top+1:returny_start-1, armx_left) ...
                                      grid(y_top+1:returny_start-1, armx_right)]')];
    
    % Plot
    plot(linear_grid)
    xline(x_ticks, LineStyle="--")
    ylabel(name)
    xlim([0,length(linear_grid)])
    xticks(x_ticks)
    xticklabels(x_labels)
    box off
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
