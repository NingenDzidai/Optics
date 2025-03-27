clc;
clear;
close all;

N = 512; % Размер сетки
step = 0.01; % Шаг сетки
x_max = step * (N / 2);
[x, y] = meshgrid(-x_max:step:x_max-step, -x_max:step:x_max-step);

Nf_values = 1:7; % Количество зон Френеля
lambda = 0.5; % Длина волны (умов.)
z = 1; % Расстояние до экрана наблюдения (умов.)

figure;
for k = 1:length(Nf_values)
    Nf = Nf_values(k);
    R = sqrt(Nf * lambda * z); % Радиус зоны Френеля

    % Создаем отверстие с Nf зонами Френеля
    aperture = double(x.^2 + y.^2 <= R^2);

    % Вычисляем дифракционную картину Френеля
    fresnel_factor = exp(1i * pi * (x.^2 + y.^2) / (lambda * z));
    U = aperture .* fresnel_factor;
    I = abs(fftshift(fft2(U))).^2;
    I = I / max(I(:));

    % Отображаем результат
    subplot(2, 4, k);
    imagesc(log(1 + I));
    colormap hot;
    colorbar;
    title(['N_f = ', num2str(Nf)]);

    % Добавляем сечения
    mid = floor(N/2);
    figure(2);
    subplot(2, 4, k);
    plot(x(mid, :), I(mid, :), 'r');
    hold on;
    plot(y(:, mid), I(:, mid), 'b');
    hold off;
    xlabel('Координата');
    ylabel('Интенсивность');
    title(['Сечение для N_f = ', num2str(Nf)]);
    legend('По x', 'По y');
end

sgtitle('Влияние количества зон Френеля на дифракционную картину');

