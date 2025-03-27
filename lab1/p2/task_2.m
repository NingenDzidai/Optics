N = 512;
step = sqrt(1/N);
x_max = step*(N/2);
[x,y] = meshgrid(-x_max:step:x_max-step, -x_max:step:x_max-step);
func = zeros(N,N);
func(257+22,:) = 1;
func(257-22,:) = 1;
Y=fftshift(func);
Y=ifft2(Y);
Y=fftshift(Y);
Y=Y*N;
Y = abs(Y).^2;
subplot(2,3,1);
pcolor(x,y, func);
colormap(gray);
axis equal;
axis([-x_max x_max -x_max x_max]);
shading interp;
title("Функция");
subplot(2,3,2:3);
plot(x(N/2+1,:),func(N/2+1,:), y(:,N/2+1),func(:,N/2+1),'r','LineWidth', 1.3);
grid on
legend("по X", "по Y")
title("Сечение функции");
subplot(2,3,4);
pcolor(x,y, abs(Y));
colormap(gray);
axis equal;
axis([-x_max x_max -x_max x_max]);
shading interp;
title("Дифракция Фраунгофера");
subplot(2,3,5:6);
plot(x(N/2+1,:),abs(Y(N/2+1,:)), y(:,N/2+1),abs(Y(:,N/2+1)),'r', 'LineWidth', 1.3);
legend("по X", "по Y")
grid on
title("Сечение дифракции Фраунгофера");
