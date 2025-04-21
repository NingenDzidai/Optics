%% Optical System Simulation with Elliptical Pupil
% Author: [Your Name]
% Date: [Current Date]
% Description: Simulates optical system with elliptical pupil and rotated input image

function elliptical_pupil_simulation()
    %% Initialization and Configuration
    clc; clear; close all;
    pkg load image;

    % Simulation parameters
    config = struct(...
        'N', 512, ...              % Grid size
        'A', 0.5, ...              % Aperture
        'lambda', 0.5, ...         % Wavelength [μm]
        'D_zr', 20, ...            % Pupil diameter
        'C', 0.7, ...              % Aberration coefficient
        'rotation_angle', 180, ... % Image rotation angle [deg]
        'pupil_a', 2, ...          % Horizontal semi-axis
        'pupil_b', 0.1, ...        % Vertical semi-axis
        'image_path', "G:/Учеба/University/10_semester/оптические_изображения/лабы/lab4/murloc.jpg"...
    );

    %% Setup Environment
    try
        % Load and prepare image
        [item, config] = load_and_prepare_image(config);

        % Initialize simulation grids and pupil
        [grids, pupil] = initialize_simulation_parameters(config);

        %% Optical Processing
        % Coherent imaging
        coherent_result = process_coherent_imaging(item, pupil, grids);

        % Incoherent imaging
        incoherent_result = process_incoherent_imaging(item, pupil, grids);

        % Calculate PSF and MTF
        [psf, mtf] = calculate_psf_mtf(pupil, grids);

        %% Visualization
        visualize_results(item, coherent_result, incoherent_result, psf, mtf, pupil, grids, config);

    catch ME
        handle_error(ME);
    end
end

%% Helper Functions
function [item, config] = load_and_prepare_image(config)
    % Load, rotate, and preprocess input image
    try
        % Read image
        pict = imread(config.image_path);

        % Rotate image
        pict = imrotate(pict, config.rotation_angle, 'bilinear', 'crop');

        % Convert to grayscale and resize
        if size(pict, 3) == 3
            pict = rgb2gray(pict);
        end
        item = double(imresize(pict, [config.N, config.N]));

    catch
        error('Image processing failed. Verify file path and format.');
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

    % Create elliptical pupil
    rho_px = rho .* sin(phi);
    rho_py = rho .* cos(phi);
    pupil = ((rho_px/config.pupil_a).^2 + (rho_py/config.pupil_b).^2) < 1;

    % Package all grids
    grids = struct(...
        'im_x', im_x, 'im_y', im_y, ...
        'p_x', p_x, 'p_y', p_y, ...
        'step_zr', step_zr, 'step_it', step_it, 'step_im', step_im, ...
        'rho_px', rho_px, 'rho_py', rho_py);
end

function result = process_coherent_imaging(item, pupil, grids)
    % Coherent imaging processing
    fft_item = (1/grids.N) * fftshift(fft2(fftshift(item)));
    result_complex = grids.N * fftshift(ifft2(fftshift(fft_item .* pupil)));
    result = abs(result_complex).^2;
end

function result = process_incoherent_imaging(item, pupil, grids)
    % Incoherent imaging processing
    fft_intens = (1/grids.N) * fftshift(fft2(fftshift(item)));
    psf_complex = fftshift(ifft2(fftshift(pupil)));
    fft_psf = 0.25*grids.N * fftshift(fft2(fftshift(abs(psf_complex).^2)));
    result_spectrum = fft_intens .* fft_psf;
    result_complex = grids.N * fftshift(ifft2(fftshift(result_spectrum)));
    result = abs(result_complex);
end

function [psf, mtf] = calculate_psf_mtf(pupil, grids)
    % Calculate PSF and MTF
    psf_complex = (grids.step_zr/grids.step_it) * (fftshift(ifft2(fftshift(pupil))) * grids.N;
    psf = (abs(psf_complex).^2) / (pi^2);

    mtf_complex = (grids.step_it/grids.step_zr) * (fftshift(fft2(fftshift(psf))) / grids.N;
    mtf = abs(mtf_complex * pi);
end

function visualize_results(item, coherent_result, incoherent_result, psf, mtf, pupil, grids, config)
    % Calculate display limits
    x_max = grids.step_im * grids.N/2;
    p_max = grids.step_zr * grids.N/2;
    center_idx = grids.N/2 + 1;

    %% Figure 1: PSF and MTF slices
    figure(1);
    set(gcf, 'Position', [100, 100, 1200, 500]);

    % PSF slice
    subplot(1,2,1);
    plot(grids.im_x(center_idx,:), psf(center_idx,:), 'b', 'LineWidth', 2);
    grid on;
    xlabel("x', μm");
    ylabel('Intensity');
    title('Point Spread Function (PSF)');
    xlim([-x_max/4, x_max/4]);

    % MTF slice
    subplot(1,2,2);
    plot(grids.p_x(center_idx,:), mtf(center_idx,:), 'b', 'LineWidth', 2);
    grid on;
    xlabel('Spatial frequency, cycles/μm');
    ylabel('Modulation');
    title('Modulation Transfer Function (MTF)');
    xlim([0, p_max/3]);

    %% Figure 2: 2D PSF
    figure(2);
    imagesc(grids.im_x(1,:), grids.im_y(:,1), psf);
    axis equal tight;
    shading interp;
    colormap('hot');
    colorbar;
    title('2D Point Spread Function');
    xlabel('x, μm');
    ylabel('y, μm');

    %% Figure 3: Original and degraded images
    figure(3);
    set(gcf, 'Position', [100, 100, 1200, 500]);

    % Original image
    subplot(1,2,1);
    imagesc(item);
    colormap('gray');
    axis equal tight;
    title('Original Image');

    % Degraded image
    subplot(1,2,2);
    imagesc(grids.im_x(1,:), grids.im_y(:,1), incoherent_result);
    colormap('gray');
    axis equal tight;
    title('Degraded Image (Incoherent)');
    xlabel('x, μm');
    ylabel('y, μm');

    %% Figure 4: Pupil shape
    figure(4);
    imagesc(grids.p_x(1,:), grids.p_y(:,1), pupil);
    axis equal tight;
    colormap('gray');
    caxis([0 1]);
    colorbar;
    title('Elliptical Pupil Shape');
    xlabel('x, pupil units');
    ylabel('y, pupil units');
end

function handle_error(ME)
    fprintf(2, '\n!!! Simulation Error !!!\n');
    fprintf(2, 'Message: %s\n', ME.message);
    fprintf(2, 'In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end

%% Run the simulation
elliptical_pupil_simulation();
