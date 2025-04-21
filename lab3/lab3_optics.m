% Инициализация
clc; clear; close all;
%%
% Загрузка и обработка изображения
try
  image = imread('rofl.jpg');
catch
  error('Could not load image. Check file path.');
end
% Определение размеров изображения
[height, width, ~] = size(image);
% Определение минимального размера стороны
minSize = min(height, width);
% Обрезка изображения до квадратного размера
croppedImage = imcrop(image, [0, 0, minSize-1, minSize-1]);
% Изменение размера изображения до [512, 512]
resizedImage = imresize(croppedImage, [512, 512]);
item = im2gray(resizedImage);
% Параметры моделирования
N = 512;
A = 0.5;
lambda = 0.5;
D_zr = 20; % 4
step_zr = D_zr / N;
step_it = 1/(N * step_zr);
step_im = step_it * lambda / A;
% Создание координатных сеток
[im_axis_X, im_axis_Y] = meshgrid(-(N/2)*step_im:step_im:(N/2-1)*step_im, ...
-(N/2)*step_im:step_im:(N/2-1)*step_im);
zrachok = zeros(N);
[p_x, p_y] = meshgrid(-(N/2)*step_zr:step_zr:(N/2-1)*step_zr, ...
-(N/2)*step_zr:step_zr:(N/2-1)*step_zr);
% Вычисление угла
fi = zeros(N);
for i = 1:N
  for j = 1:N
    if (p_y(i,j) == 0)
      if (p_x(i,j) >= 0)
        fi(i,j) = pi / 2;
      else
        fi(i,j) = -pi / 2;
      end
    else
      if (p_y(i,j) >= 0)
        fi(i,j) = atan(p_x(i,j) / p_y(i,j));
      else
        fi(i,j) = pi + atan(p_x(i,j) / p_y(i,j));
      end
    end
  end
end
% Расчет параметров
ro = sqrt(p_x.*p_x + p_y.*p_y);
ro_p_x = ro.*sin(fi);
ro_p_y = ro.*cos(fi);
% Задание аберрации
C = 0.455;
aberr = exp(2*pi*1i*C*(3*ro.^3 - 2*ro).*cos(fi));
% Формирование зрачка
for i = 1:N
  for j = 1:N
    if sqrt(ro_p_x(i,j)^2 + ro_p_y(i,j)^2) < 1
      zrachok(i,j) = 1;
    end
  end
end
zrachok = zrachok.*aberr;
% Когерентное изображение
fft_item = 1/N * (fftshift(fft2(fftshift(item))));
res = N * (fftshift(ifft2(fftshift(fft_item.*zrachok))));
res = abs(res).^2;
res = flipud(res);
% Некогерентное изображение
fft_intens = 1/N * (fftshift(fft2(fftshift(abs(item)))));
fft_func_rasp = 0.25*N * fftshift(fft2(fftshift(abs(fftshift(ifft2(fftshift(zrachok)))).^2)));
func_rasp_img = fft_intens .* fft_func_rasp;
intens_rasp_img = N * (fftshift(ifft2(fftshift(func_rasp_img))));
intens_rasp_img = flipud(intens_rasp_img);
% Подготовка к визуализации
[x, y] = meshgrid(-(N/2)*step_im:step_im:(N/2-1)*step_im, ...
-(N/2)*step_im:step_im:(N/2-1)*step_im);
n_max = step_it*N/2;
x_max = step_im*N/2;
p_max = step_zr*N/2;
% Расчет ФРТ и ФПМ
FRT_ = (step_zr/step_it)*(fftshift(ifft2(fftshift(zrachok)))*N);
FRT_abs = (abs(FRT_).*abs(FRT_))/(pi^2);
D = (step_it/step_zr)*(fftshift(fft2(fftshift(FRT_abs)))/N);
D_norm = D*pi;
D_abs = abs(D_norm);
% Вывод числа Штреля
[max_value, max_index] = max(FRT_abs(:, N/2+1));
fprintf('Число Штреля: %f\n', max_value);
% Визуализация результатов
figure('Position', [500, 150, 1000, 800]);
sgtitle(sprintf('Коэффициент разложения: C = %.3f\nЧисло Штреля: %.3f', C, max_value), 'FontSize', 16)
% Срез ФРТ
subplot(2,2,1)
plot(x(N/2+1,:), FRT_abs(:, N/2+1), 'Color', 'r', 'LineWidth', 1.3)
xlim([-x_max/4, x_max/4]) % Оставляем пределы по оси X
grid on
title('Срез ФРТ')
% ФПМ
subplot(2,2,2)
plot(p_x(N/2+1,:), D_abs(N/2+1,:), 'Color', 'r', 'LineWidth', 1.3)
xlim([0, p_max/3]) % Оставляем пределы по оси X
xticks(0:0.2:p_max/3)
grid on
title('ФПМ')
% Полноценная ФРТ
subplot(2,2,3)
pcolor(im_axis_X, im_axis_Y, FRT_abs)
colormap('gray')
axis equal tight
shading interp
xlim([-x_max, x_max])
title('Функция рассеяния точки')
% Некогерентное изображение
subplot(2,2,4)
pcolor(im_axis_X, im_axis_Y, intens_rasp_img)
colormap('gray')
axis equal tight
shading interp
xlim([-x_max, x_max])
title('Некогерентное изображение')
