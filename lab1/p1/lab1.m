clc; clear; close all;

% Определяем параметры
Nx = 1000;  % Количество точек по x
Ny = 1000;  % Количество точек по y
Lx = 4;    % Размер области в пространстве отверстия (по x)
Ly = 2;    % Размер области в пространстве отверстия (по y)

x = linspace(-Lx/2, Lx/2, Nx);
y = linspace(-Ly/2, Ly/2, Ny);
[X, Y] = meshgrid(x, y);

% Задаем отверстие rect(2x) * rect(y)
aperture = double(abs(X) <= 1 & abs(Y) <= 0.5);

figure;
imagesc(x, y, aperture);
colormap gray;
axis equal;
title('Апертура (прямоугольное отверстие)');


% Вычисляем двумерное Фурье-преобразование (FFT)
I = abs(fftshift(fft2(aperture))).^2;

% Определяем координаты в пространстве Фурье
fx = linspace(-Nx/(2*Lx), Nx/(2*Lx), Nx);
fy = linspace(-Ny/(2*Ly), Ny/(2*Ly), Ny);
[FX, FY] = meshgrid(fx, fy);

% Нормализация интенсивности
I = I / max(I(:));

% Визуализация
figure;
imagesc(fx, fy, I);
colormap jet;
colorbar;
xlabel('fx');
ylabel('fy');
title('Дифракционная картина Фраунгофера (интенсивность)');

