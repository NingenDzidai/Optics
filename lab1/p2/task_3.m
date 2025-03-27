N = 512;
r=0.051;
x_max=0.256;
step = x_max*2/N;
lambda=0.5e-6;
z=10;
n=1;
func = zeros(N,N);
[x,y] = meshgrid(-x_max:step:x_max-step, -x_max:step:x_max-step);
k=2*pi/lambda;
for i=1:1:N
  for j=1:1:N
    if (power(x(i,j),2) + power(y(i,j),2)) < power(r, 2);
      func(i,j) = 1;
    end
  end
end
h = zeros(N, N);
for i=1:1:N
  for j=1:1:N
    h(i,j) = 1/(z*lambda*1i) * exp(1i*k*n*z) * exp(1i*k*n/(2*z)*(x(i,j)^2+y(i,j)^2));
  end
end
sv = conv2(func, h);
sv = abs(sv(257:768, 257:768));
subplot(2,3,1);
pcolor(x,y, func);
colormap(gray);
axis equal;
axis([-x_max x_max -x_max x_max]);
shading interp;
title("Функция");
subplot(2,3,2:3);
plot(x(N/2+1,:),func(N/2+1,:), y(:,N/2+1),func(:,N/2+1),'r', 'LineWidth', 1.3);
title("Сечение функции");
legend("по X", "по Y");
subplot(2,3,4);
pcolor(x,y, abs(sv));
colormap(gray);
axis equal;
axis([-x_max x_max -x_max x_max]);
shading interp;
title("Дифракция Френеля");
subplot(2,3,5:6);
plot(x(N/2+1,:),abs(sv(N/2+1,:)), y(:,N/2+1),abs(sv(:,N/2+1)),'r','LineWidth', 1.3);
title("Сечение дифракции Френеля");
legend("по X", "по Y");
fresnel_number = (r^2)/(lambda*z);
disp(['Fresnel number: ' num2str(fresnel_number)]);
