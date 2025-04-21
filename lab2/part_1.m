% Часть 1. Формирование изображения предмета-решётки
clear;
clc;
pkg load image;  % Для некоторых версий Octave может потребоваться

% Параметры системы
N = 512;         % Размер массива
lambda = 0.5;    % Длина волны, мкм
Dzr = 20;        % Охват зрачка
A = 0.5;         % Апертура

% Вычисление производных параметров
Rzr = N / Dzr;               % Радиус зрачка
dp = Dzr / N;                % Шаг по зрачку, к. е.
dn = 1 / (N * dp);           % Шаг по предмету, к. е.
dx = dn * (lambda / A);      % Шаг по изображению, мкм
n_max = dn * (N/2);          % Максимальное значение координаты предмета
p_max = dp * (N/2);          % Максимальное значение координаты зрачка
x_max = dx * (N/2);          % Максимальное значение координаты изображения

% Создание массива координат
coords = linspace(-n_max, n_max-dn, N);
[nx, ny] = meshgrid(coords, coords);

coords_p = linspace(-p_max, p_max-dp, N);
[px, py] = meshgrid(coords_p, coords_p);

coords_x = linspace(-x_max, x_max-dx, N);
[x, y] = meshgrid(coords_x, coords_x);

% Функция зрачка
f = (px.^2 + py.^2) < 1;  % Rzr = 1, к. е.

% Обработка решетки с заданным числом линий (count)
for count = 1:2:17
    % Генерация решетки
    w = round(N/(count*2));  % Ширина линии
    grid = zeros(N,N);

    % Создаем одну линию решетки
    line_pattern = [zeros(1,w), ones(1,w)];
    full_pattern = repmat(line_pattern, 1, ceil(N/(2*w)));
    full_pattern = full_pattern(1:N);

    % Симметризуем
    full_pattern(1+N/2:end) = fliplr(full_pattern(1:N/2));

    % Создаем 2D решетку
    item = repmat(full_pattern, N, 1);

    % Визуализация предмета
    subplot(2, 3, 1);
    imagesc(x(1,:), y(:,1), item);
    colormap(gray);
    axis equal tight;
    xlabel('x, мкм');
    ylabel('y, мкм');
    title('Предмет');

    % Когерентное освещение
    fft_item = fft2(item) / N;
    filtered_image = fft_item .* f;
    ifft_image = ifft2(filtered_image) * N;
    intensity_image = abs(ifft_image).^2;

    subplot(2, 3, 2);
    imagesc(x(1,:), y(:,1), intensity_image);
    colormap(gray);
    axis equal tight;
    xlabel("x', мкм");
    ylabel("y', мкм");
    title('Когерентное освещение');

    % Некогерентное освещение
    object_intensity = abs(item).^2;
    fft_object = fft2(object_intensity);

    PSF = abs(ifft2(f)).^2;
    PSF = PSF / sum(PSF(:));  % Нормировка
    OTF = fft2(PSF);

    fft_image = fft_object .* OTF;
    ifft_intens_image = real(ifft2(fft_image));
    ifft_intens_image = ifft_intens_image / max(ifft_intens_image(:));

    subplot(2, 3, 3);
    imagesc(x(1,:), y(:,1), ifft_intens_image);
    colormap(gray);
    axis equal tight;
    xlabel("x', мкм");
    ylabel("y', мкм");
    title('Некогерентное освещение');

    % Визуализация сечений
    subplot(2, 1, 2);
    plot(x(N/2, :), item(N/2, :), 'k-', 'LineWidth', 2, ...
         x(N/2, :), intensity_image(N/2, :), 'r-', 'LineWidth', 2, ...
         x(N/2, :), ifft_intens_image(N/2, :), 'b-', 'LineWidth', 2);
    xlim([-x_max, x_max]);
    grid on;
    legend('Предмет', 'Когерентное освещение', 'Некогерентное освещение');
    xlabel("x', мкм");
    title('Сечение');

    % Общий заголовок
    sgtitle(['Количество линий: ' num2str(count)], 'FontSize', 14);
    pause(1);
end
