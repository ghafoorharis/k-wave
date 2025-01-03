function [scan_line_example,t0,r,scan_lines_fund,scan_lines_harm] = process_results(kgrid, input_signal, ...
                                medium, scan_lines, ...
                                tone_burst_freq, ...
                                c0,number_scan_lines)
    scan_line_win = getWin(kgrid.Nt * 2, 'Tukey', 'Param', 0.05).';
    scan_line_win = [zeros(1, length(input_signal) * 2), scan_line_win(1:end/2 - length(input_signal) * 2)];
    scan_lines = bsxfun(@times, scan_line_win, scan_lines);
    scan_line_example(1, :) = scan_lines(end/2, :);
    
    % Time Gain Compensation
    t0 = length(input_signal) * kgrid.dt / 2;
    r = c0 * ( (1:length(kgrid.t_array)) * kgrid.dt - t0 ) / 2;
    tgc_alpha_db_cm = medium.alpha_coeff * (1.5e6 * 1e-6)^medium.alpha_power;
    tgc = exp(tgc_alpha_db_cm / 8.686 * 100 * 2 * r);
    scan_lines = bsxfun(@times, tgc, scan_lines);
    scan_line_example(2, :) = scan_lines(end/2, :);
    
    % Frequency Filtering
    scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
    scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);
    scan_line_example(3, :) = scan_lines_fund(end/2, :);
    
    % Envelope Detection
    scan_lines_fund = envelopeDetection(scan_lines_fund);
    scan_lines_harm = envelopeDetection(scan_lines_harm);
    scan_line_example(4, :) = scan_lines_fund(end/2, :);
    
    % Log Compression
    compression_ratio = 30;
    scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
    scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);
    scan_line_example(5, :) = scan_lines_fund(end/2, :);
    scan_line_example(6, :) = scan_lines_harm(end/2, :);
    % upsample the image using linear interpolation
    scale_factor = 2;
    scan_lines_fund = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
    scan_lines_harm = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_harm, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');

end
