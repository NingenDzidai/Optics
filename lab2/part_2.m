% Часть 2. Формирование изображения произвольного изображения
% Полностью рабочий код для Octave

clear; clc; close all;
pkg load image;

%% Параметры системы
N = 512;                    % Размер массива [пиксели]
lambda = 0.5;               % Длина волны [мкм]
Dzr = 20;                   % Охват зрачка [к.е.]
A = 0.5;                    % Апертура [к.е.]
Rzr = N / Dzr;              % Радиус зрачка [пиксели]
dp = Dzr / N;               % Шаг по зрачку [к.е./пиксель]
dn = 1 / (N * dp);          % Шаг по предмету [к.е.]
dx = dn * (lambda / A);     % Шаг по изображению [мкм]

%% Создание координатных сеток
n_values = linspace(-dn*N/2, dn*N/2-dn, N);
p_values = linspace(-dp*N/2, dp*N/2-dp, N);
x_values = linspace(-dx*N/2, dx*N/2-dx, N);

[nx, ny] = meshgrid(n_values, n_values);
[px, py] = meshgrid(p_values, p_values);
[x, y] = meshgrid(x_values, x_values);

%% Функция зрачка
r = sqrt(px.^2 + py.^2);
f = double(r <= 1); % Круглая апертура

%% Загрузка и подготовка изображения
[fn, pn] = uigetfile({'*.jpg;*.png;*.bmp;*.tif', 'Image Files'}, 'Select an image');
if isequal(fn, 0)
    item = double(zeros(N));
    item(N/4:3*N/4, N/4:3*N/4) = 1;
    item(N/3:2*N/3, N/3:2*N/3) = 0;
else
    img = im2double((imread(fullfile(pn, fn))));
    item = imresize(img, [N, N]);
end

%% Когерентное освещение
fft_item = fftshift(fft2(ifftshift(item)));
filtered_image = fft_item .* f;
ifft_image = fftshift(ifft2(ifftshift(filtered_image)));
intensity_image = abs(ifft_image).^2;

%% Некогерентное освещение
object_intensity = abs(item).^2;
fft_object = fftshift(fft2(ifftshift(object_intensity)));

% Расчет PSF и OTF
ifft_scattering = (dp/dn) * N * fftshift(ifft2(ifftshift(f)));
PSF = abs(ifft_scattering).^2 / (pi^2);
OTF = (dn/dp) * pi * fftshift(fft2(ifftshift(PSF))) / N;

fft_image = fft_object .* OTF;
ifft_intens_image = real(fftshift(ifft2(ifftshift(fft_image))));
ifft_intens_image = ifft_intens_image / max(ifft_intens_image(:));

%% Визуализация результатов
figure('Name', 'Анализ формирования изображения', 'Position', [100 100 1200 800]);

% 1. Исходное изображение
subplot(2,3,1);
imagesc(x_values, x_values, item);
colormap(gray); axis equal tight; colorbar;
title('Исходное изображение');
xlabel('x, мкм'); ylabel('y, мкм');

% 2. Когерентное изображение
subplot(2,3,2);
imagesc(x_values, x_values, intensity_image);
colormap(gray); axis equal tight; colorbar;
title('Когерентное изображение');
xlabel('x, мкм'); ylabel('y, мкм');

% 3. Некогерентное изображение
subplot(2,3,3);
imagesc(x_values, x_values, ifft_intens_image);
colormap(gray); axis equal tight; colorbar;
title('Некогерентное изображение');
xlabel('x, мкм'); ylabel('y, мкм');

% 4. Фаза когерентного изображения
subplot(2,3,4);
imagesc(x_values, x_values, angle(ifft_image));
colormap(jet); axis equal tight; colorbar;
title('Фаза когерентного изображения');
xlabel('x, мкм'); ylabel('y, мкм');

% 5. Сравнение сечений
subplot(2,3,[5,6]);
plot(x_values, item(N/2,:), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Исходное');
hold on;
plot(x_values, intensity_image(N/2,:)/max(intensity_image(:)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Когерентное');
plot(x_values, ifft_intens_image(N/2,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Некогерентное');
hold off;
grid on; xlim([-x_max x_max]); ylim([0 1.1]);
legend('Location', 'northeast');
title('Сравнение сечений (нормированные)');
xlabel('x, мкм'); ylabel('Нормированная интенсивность');

%% Дополнительное окно с характеристиками системы
figure('Name', 'Характеристики оптической системы', 'Position', [200 200 1000 600]);

% 1. Функция зрачка
subplot(2,3,1);
imagesc(p_values, p_values, f);
colormap(gray); axis equal tight; colorbar;
title('Функция зрачка');
xlabel('px, к.е.'); ylabel('py, к.е.');

% 2. ФРТ (PSF)
subplot(2,3,2);
imagesc(x_values, x_values, PSF/max(PSF(:)));
colormap(gray); axis equal tight; colorbar;
title('ФРТ (нормированная)');
xlabel('x, мкм'); ylabel('y, мкм');

% 3. Сечение ФРТ
subplot(2,3,3);
plot(x_values, PSF(N/2,:)/max(PSF(:)), 'b-', 'LineWidth', 2);
xlim([-x_max/4 x_max/4]); grid on;
title('Сечение ФРТ');
xlabel('x, мкм'); ylabel('Нормированная интенсивность');

% 4. Модуль ФПМ (OTF)
subplot(2,3,4);
plot(p_values, abs(OTF(N/2,:)), 'r-', 'LineWidth', 2);
xlim([0 p_max/2]); ylim([0 1.1]); grid on;
title('Модуль ФПМ');
xlabel('Пространственная частота, к.е.'); ylabel('Модуль передачи');

% 5. Фаза ФПМ
subplot(2,3,5);
plot(p_values, angle(OTF(N/2,:)), 'g-', 'LineWidth', 2);
xlim([0 p_max/2]); grid on;
title('Фаза ФПМ');
xlabel('Пространственная частота, к.е.'); ylabel('Фаза [рад]');

% 6. 3D ФРТ (исправленная версия)
subplot(2,3,6);
step = 8; % Уменьшаем количество точек для 3D графика
idx = 1:step:N;
[x_plot, y_plot] = meshgrid(x_values(idx), y_values(idx));
surf(x_plot, y_plot, PSF(idx,idx)/max(PSF(:)));
shading interp; colormap(jet);
title('3D вид ФРТ');
xlabel('x, мкм'); ylabel('y, мкм'); zlabel('Интенсивность');
view(30, 45);

%% Сохранение результатов
save_results = questdlg('Сохранить результаты?', 'Сохранение', 'Да', 'Нет', 'Да');
if strcmp(save_results, 'Да')
    folder_name = uigetdir('', 'Выберите папку для сохранения');
    if folder_name ~= 0
        saveas(gcf, fullfile(folder_name, 'optical_system_results.png'));
        f = findobj('Name', 'Характеристики оптической системы');
        if ~isempty(f)
            saveas(f, fullfile(folder_name, 'system_characteristics.png'));
        end
        save(fullfile(folder_name, 'simulation_data.mat'), ...
             'item', 'intensity_image', 'ifft_intens_image', 'PSF', 'OTF');
        disp('Результаты успешно сохранены');
    end
end
