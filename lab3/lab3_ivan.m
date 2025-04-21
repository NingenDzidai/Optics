%% Optical System Simulation with 3rd Order Astigmatism Aberration
% Author: [Your Name]
% Date: [Current Date]
% Description: Simulates the effects of 3rd order Zernike astigmatism on
% PSF, MTF, and image quality under incoherent illumination.

function main()
    %% Initialization and Configuration
    clc; clear; close all;

    % Simulation parameters
    config = struct(...
        'N', 512, ...              % Grid size
        'A', 0.5, ...              % Aperture
        'lambda', 0.5, ...         % Wavelength [μm]
        'D_zr', 20, ...            % Pupil diameter
        'save_dir', 'sim_results', % Output directory
        'image_file', 'rofl.jpg'...% Input image
    );

    %% Setup Environment
    try
        [item, config] = load_and_prepare_image(config);
        [grids, pupil] = initialize_simulation_parameters(config);

        % Find optimal aberration coefficients
        [C_values, metrics] = find_optimal_coefficients(pupil, grids, config);

        % Visualize all cases
        visualize_all_cases(C_values, pupil, item, grids, config);

    catch ME
        handle_error(ME);
    end
end

%% Helper Functions
function [item, config] = load_and_prepare_image(config)
    % Load and preprocess input image
    try
        img = imread(config.image_file);
        if size(img, 3) == 3
            img = rgb2gray(img);
        end
        item = double(imresize(img, [config.N, config.N]));

        % Create output directory
        if ~exist(config.save_dir, 'dir')
            mkdir(config.save_dir);
            fprintf('Created output directory: %s\n', config.save_dir);
        end
    catch
        error('Image loading failed. Verify file existence and format.');
    end
end

function [grids, pupil] = initialize_simulation_parameters(config)
    % Calculate sampling steps
    step_zr = config.D_zr / config.N;
    step_it = 1/(config.N * step_zr);
    step_im = step_it * config.lambda / config.A;

    % Create coordinate grids
    [im_x, im_y] = meshgrid(...
        linspace(-config.N/2, config.N/2-1, config.N) * step_im, ...
        linspace(-config.N/2, config.N/2-1, config.N) * step_im);

    [p_x, p_y] = meshgrid(...
        linspace(-config.N/2, config.N/2-1, config.N) * step_zr, ...
        linspace(-config.N/2, config.N/2-1, config.N) * step_zr);

    % Polar coordinates
    [phi, rho] = cart2pol(p_x, p_y);
    rho_norm = rho / (config.D_zr/2); % Normalized radius

    % Create pupil function
    pupil = double(rho_norm <= 1);

    % Package all grids
    grids = struct(...
        'im_x', im_x, 'im_y', im_y, ...
        'p_x', p_x, 'p_y', p_y, ...
        'rho_norm', rho_norm, 'phi', phi, ...
        'step_zr', step_zr, 'step_it', step_it, 'step_im', step_im);
end

function [C_values, metrics] = find_optimal_coefficients(pupil, grids, config)
    % Target metrics to optimize
    targets = struct(...
        'Strehl_08', 0.8, ...
        'Strehl_05', 0.5, ...
        'Contrast_s10', 0.2, ...
        'Contrast_s08', 0.2);

    % Initialize results
    C_values = struct();
    metrics = struct();

    fprintf('\n=== Optimizing Aberration Coefficients ===\n');

    % Find coefficients for each target
    C_values.Strehl_08 = optimize_coefficient(...
        'Strehl', targets.Strehl_08, pupil, grids, config);
    C_values.Strehl_05 = optimize_coefficient(...
        'Strehl', targets.Strehl_05, pupil, grids, config);
    C_values.Contrast_s10 = optimize_coefficient(...
        'Contrast_s10', targets.Contrast_s10, pupil, grids, config);
    C_values.Contrast_s08 = optimize_coefficient(...
        'Contrast_s08', targets.Contrast_s08, pupil, grids, config);

    % Display results
    fprintf('\n=== Optimal Coefficients Found ===\n');
    fprintf('Strehl ≈ 0.8: C = %.4f\n', C_values.Strehl_08);
    fprintf('Strehl ≈ 0.5: C = %.4f\n', C_values.Strehl_05);
    fprintf('Contrast 0.2 @ s≈1.0: C = %.4f\n', C_values.Contrast_s10);
    fprintf('Contrast 0.2 @ s≈0.8: C = %.4f\n', C_values.Contrast_s08);
end

function C_opt = optimize_coefficient(target_type, target_value, pupil, grids, config)
    % Optimization parameters
    search_params = struct(...
        'C_min', 0.01, ...
        'C_max', 2.0, ...
        'coarse_steps', 100, ...
        'fine_steps', 20, ...
        'tolerance', 0.05);

    fprintf('\nOptimizing for %s = %.2f\n', target_type, target_value);

    % Coarse search
    C_range = linspace(search_params.C_min, search_params.C_max, search_params.coarse_steps);
    [C_opt, min_diff] = evaluate_candidates(C_range, target_type, target_value, pupil, grids);

    % Fine search
    C_range = linspace(...
        max(search_params.C_min, C_opt - search_params.tolerance), ...
        min(search_params.C_max, C_opt + search_params.tolerance), ...
        search_params.fine_steps);
    [C_opt, ~] = evaluate_candidates(C_range, target_type, target_value, pupil, grids);

    fprintf('Optimal C = %.4f\n', C_opt);
end

function [C_opt, min_diff] = evaluate_candidates(C_range, target_type, target_value, pupil, grids)
    min_diff = Inf;
    C_opt = 0;

    for C = C_range
        metrics = calculate_metrics(C, pupil, grids);

        switch target_type
            case 'Strehl'
                current_value = metrics.strehl;
            case 'Contrast_s10'
                current_value = metrics.contrast_s10;
            case 'Contrast_s08'
                current_value = metrics.contrast_s08;
        end

        diff = abs(current_value - target_value);

        if diff < min_diff
            min_diff = diff;
            C_opt = C;
        end
    end
end

function metrics = calculate_metrics(C, pupil, grids)
    % Apply aberration
    W = C * grids.rho_norm.^2 .* cos(2*grids.phi);
    pupil_aberr = pupil .* exp(1i * 2*pi * W);

    % Calculate PSF
    PSF = abs(fftshift(ifft2(ifftshift(pupil_aberr)))).^2;
    PSF = PSF / sum(PSF(:)); % Energy normalization

    % Calculate MTF
    MTF = abs(fftshift(fft2(ifftshift(PSF))));
    MTF = MTF / max(MTF(:)); % Normalize to DC

    % Calculate metrics
    center = grids.N/2 + 1;
    metrics.strehl = PSF(center, center);

    % Frequency indices
    freq_08 = round(0.8 * center/2) + center;
    freq_10 = round(1.0 * center/2) + center;
    freq_08 = min(freq_08, grids.N);
    freq_10 = min(freq_10, grids.N);

    metrics.contrast_s08 = MTF(center, freq_08);
    metrics.contrast_s10 = MTF(center, freq_10);
end

function visualize_all_cases(C_values, pupil, item, grids, config)
    fprintf('\n=== Generating Visualization ===\n');

    % Case 1: No aberration
    visualize_case(0, pupil, item, grids, config, 'No Aberration');

    % Other cases
    visualize_case(C_values.Strehl_08, pupil, item, grids, config, 'Strehl ≈ 0.8');
    visualize_case(C_values.Strehl_05, pupil, item, grids, config, 'Strehl ≈ 0.5');
    visualize_case(C_values.Contrast_s10, pupil, item, grids, config, 'Contrast 0.2 @ s≈1.0');
    visualize_case(C_values.Contrast_s08, pupil, item, grids, config, 'Contrast 0.2 @ s≈0.8');

    fprintf('All visualizations saved to: %s\n', config.save_dir);
end

function visualize_case(C, pupil, item, grids, config, case_name)
    % Calculate metrics and aberrated fields
    metrics = calculate_metrics(C, pupil, grids);

    % Apply aberration
    W = C * grids.rho_norm.^2 .* cos(2*grids.phi);
    pupil_aberr = pupil .* exp(1i * 2*pi * W);

    % Calculate PSF and MTF
    PSF = abs(fftshift(ifft2(ifftshift(pupil_aberr)))).^2;
    PSF = PSF / max(PSF(:));
    MTF = abs(fftshift(fft2(ifftshift(PSF))));
    MTF = MTF / max(MTF(:));

    % Calculate degraded image
    image_spectrum = fftshift(fft2(item));
    degraded_spectrum = image_spectrum .* MTF;
    degraded_image = abs(ifft2(ifftshift(degraded_spectrum)));

    % Create figure
    fig = figure('Position', [100, 100, 1200, 900], 'Color', 'w');

    % Main title
    annotation('textbox', [0.1, 0.95, 0.8, 0.05], ...
        'String', sprintf('%s (C = %.3f, Strehl = %.3f)', case_name, C, metrics.strehl), ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold');

    % PSF Slice
    subplot(2,2,1);
    plot(grids.im_x(grids.N/2+1,:), PSF(grids.N/2+1,:), 'LineWidth', 2);
    xlim([-grids.step_im*grids.N/8, grids.step_im*grids.N/8]);
    title('PSF Cross-Section');
    xlabel('Position [μm]'); ylabel('Normalized Intensity');
    grid on; box on;

    % MTF
    subplot(2,2,2);
    spatial_freq = linspace(0, 1/(2*grids.step_im), grids.N/2);
    plot(spatial_freq, MTF(grids.N/2+1, grids.N/2+1:end), 'LineWidth', 2);
    title('Modulation Transfer Function');
    xlabel('Spatial Frequency [cycles/μm]'); ylabel('Contrast');
    grid on; box on;

    % PSF 2D
    subplot(2,2,3);
    imagesc(grids.im_x(1,:), grids.im_y(:,1), PSF);
    axis image; colormap('hot'); colorbar;
    title('Point Spread Function');
    xlabel('x [μm]'); ylabel('y [μm]');

    % Degraded Image
    subplot(2,2,4);
    imagesc(degraded_image);
    axis image; colormap('gray'); colorbar;
    title('Degraded Image');

    % Save figure
    filename = regexprep(lower(case_name), {' ', '≈', '@'}, {'_', '', ''});
    saveas(fig, fullfile(config.save_dir, sprintf('%s_C%.3f.png', filename, C)));
    close(fig);
end

function handle_error(ME)
    fprintf(2, '\n!!! Simulation Error !!!\n');
    fprintf(2, 'Message: %s\n', ME.message);
    fprintf(2, 'In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);

    % Try to save partial results if possible
    try
        if exist('config', 'var') && isfield(config, 'save_dir')
            save(fullfile(config.save_dir, 'error_dump.mat'));
            fprintf(2, 'Partial results saved to error_dump.mat\n');
        end
    catch
    end
end

%% Run the simulation
main();
