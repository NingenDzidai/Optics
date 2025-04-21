pkg load image

% Initialization
clc; clear; close all;

%% Image Processing
try
    original_image = imread('rofl.jpg');
    if size(original_image, 3) == 3
        test_image = rgb2gray(imresize(original_image, [512 512]));
    else
        test_image = imresize(original_image, [512 512]);
    end
catch
    error('Could not load image. Check file path.');
end

%% Simulation Parameters
N = 512;                % Grid size
lambda = 0.5e-6;        % Wavelength [m]
D = 20e-3;              % Pupil diameter [m]
f = 0.5;                % Focal length [m]

% Spatial sampling
delta_zr = D/N;         % Pupil plane sampling [m]
delta_im = lambda*f/D;  % Image plane sampling [m]

% Coordinate grids
[im_x, im_y] = meshgrid((-N/2:N/2-1)*delta_im, (-N/2:N/2-1)*delta_im);
[p_x, p_y] = meshgrid((-N/2:N/2-1)*delta_zr, (-N/2:N/2-1)*delta_zr);

% Normalized pupil coordinates
rho = sqrt(p_x.^2 + p_y.^2)/(D/2); % Normalized to pupil radius
phi = atan2(p_y, p_x);

%% Zernike Astigmatism Aberration
C_values = [0, 0.2, 0.4]; % Aberration coefficients in wavelengths

% Initialize results
results = struct();

for i = 1:length(C_values)
    C = C_values(i);

    % Create pupil function with proper normalization
    W = C*lambda * rho.^2 .* cos(2*phi); % Wavefront aberration [m]
    pupil = double(rho <= 1) .* exp(1i * 2*pi * W/lambda);

    %% PSF Calculation (properly normalized)
    PSF_complex = fftshift(ifft2(ifftshift(pupil)));
    PSF = abs(PSF_complex).^2;
    PSF = PSF / sum(PSF(:)); % Energy normalization

    %% OTF Calculation
    OTF_complex = fftshift(fft2(ifftshift(PSF)));
    OTF = abs(OTF_complex);
    OTF = OTF / max(OTF(:)); % Normalize to DC component

    %% Image Degradation
    image_spectrum = fftshift(fft2(double(test_image)));
    degraded_spectrum = image_spectrum .* OTF;
    degraded_image = abs(ifft2(ifftshift(degraded_spectrum)));

    % Store results
    results(i).C = C;
    results(i).PSF = PSF;
    results(i).OTF = OTF;
    results(i).degraded_image = degraded_image;

    % Calculate Strehl ratio
    if i == 1
        PSF_perfect = PSF;
    end
    strehl_ratio = max(PSF(:))/max(PSF_perfect(:));
    results(i).strehl = strehl_ratio;

    fprintf('Case %d: C = %.2fλ, Strehl ratio = %.3f\n', i, C, strehl_ratio);
end

%% Visualization
for i = 1:length(results)
    figure('Position', [100, 100, 1200, 900]);
    ha = axes('Position',[0 0 1 1],'Visible','off');
    text(0.5, 0.98, sprintf('Astigmatism Aberration (Z2^2)\nC = %.2fλ, Strehl Ratio = %.3f', ...
                   results(i).C, results(i).strehl), 'FontSize', 14);

    % PSF
    subplot(2,2,1);
    imagesc(im_x(1,:)*1e6, im_y(:,1)*1e6, results(i).PSF/max(results(i).PSF(:)));
    axis image; colormap('hot'); colorbar;
    xlabel('x [μm]'); ylabel('y [μm]');
    title('Normalized Point Spread Function');
    clim([0 1]); % Ensure full dynamic range

    % PSF Slice
    subplot(2,2,2);
    plot(im_x(N/2+1,:)*1e6, results(i).PSF(N/2+1,:)/max(results(i).PSF(:)), 'b', 'LineWidth', 2);
    xlabel('Position [μm]'); ylabel('Normalized Intensity');
    title('PSF Cross-Section'); grid on; xlim([-50, 50]);
    ylim([0 1.1]);

    % OTF
    subplot(2,2,3);
    spatial_freq = (-N/2:N/2-1)/(N*delta_im);
    freq_mask = spatial_freq >= 0;
    plot(spatial_freq(freq_mask)*1e-3, results(i).OTF(N/2+1,freq_mask), 'r', 'LineWidth', 2);
    xlabel('Spatial Frequency [cycles/mm]'); ylabel('Modulation');
    title('Normalized OTF'); grid on; xlim([0 100]);
    ylim([0 1.1]);

    % Degraded Image
    subplot(2,2,4);
    imagesc(results(i).degraded_image); colormap('gray'); axis image off;
    title('Degraded Image');
end
