clear all;
close all;

%a See slide 93

%b 
x = [1,1,1,0.5];
nx = 0:3;

h = [1,0.75,0.5,0.25];
nh = 0:3;

yp = cconv(x,h,4);
subplot(2,1,1);
stem(nh,yp);
xlim([-1 7]);
title('Circular convolution result');
xlabel('n');

%c
[ny,y] = convcool(nx,x,nh,h);
subplot(2,1,2);
stem(ny,y);
xlim([-1,7]);
title('Linear convolution result');
xlabel('n');