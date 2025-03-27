clc;
clear all;
N = 512;
step = sqrt(1/N);
x_max = step*(N/2);
func = zeros(N,N);

[x,y] = meshgrid(-x_max:step:x_max-step, -x_max:step:x_max-step);

% Задаем прямоугольное отверстие rect(2x) * rect(y)
for i = 1:N
    for j = 1:N
        if (abs(x(i, j)) <= 1) && (abs(y(i, j)) <= 0.5)
            func(i, j) = 1;
        end
    end
end


Y=fftshift(func)/N;
Y=ifft2(Y);
Y=fftshift(Y);
Y=Y;
intens=abs(Y).*abs(Y);
intens=Y.*conj(Y);
subplot(2,3,1);
pcolor(x,y, func); % отрисовка 3D-карты уровней
colormap(gray); % черно-белая палитра
axis equal; % одинаковый масштаб по осям
axis([-x_max x_max -x_max x_max]); % диапазон значений по осям
shading interp; % раскраска с использование интерполяции
title("Функция");
subplot(2,3,4);
plot(x(N/2+1,:),func(N/2+1,:), y(:,N/2+1),func(:,N/2+1),'r', 'LineWidth', 1.3);
title("Сечение функции");
legend("по X", "по Y");
subplot(2,3,2);
pcolor(x,y, abs(Y));
colormap(gray);
axis equal;
axis([-x_max x_max -x_max x_max]);
shading interp;
title("Результат Фурье-преобразования");
subplot(2,3,3);
pcolor(x,y, abs(intens));
colormap(gray);
axis equal;
axis([-x_max x_max -x_max x_max]);
shading interp;
title("Распределение интенсивности");
subplot(2,3,5);
plot(x(N/2+1,:),abs(Y(N/2+1,:)), y(:,N/2+1),abs(Y(:,N/2+1)),'r','LineWidth', 1.3);
title("Сечение Фурье-преобразования");
legend("по X", "по Y");
subplot(2,3,6);
plot(x(N/2+1,:),abs(intens(N/2+1,:)), y(:,N/2+1),abs(intens(:,N/2+1)),'r','LineWidth', 1.3);
title("Сечение интенсивности");
legend("по X", "по Y");
