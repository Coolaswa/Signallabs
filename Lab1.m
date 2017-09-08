close all; clear all;

%% Assignment 1: Convolution
nx = [-3 -2 -1 0 1 2 3];
x = [1 0 1 0 1 0 1]; %start at n=-3, ends at n=3
h = [-0.5 1 -0.5 ];
y = conv(x,h);
ny = linspace(nx(1) - 1, length(y)+nx(1)-1-1,length(y));
stem(nx,x);
hold on;
stem(ny,y);
xlim([min([nx,ny])-1 max([nx,ny])+1]);
ylim([min([x,y])-0.5 max([x,y])+0.5]);
grid on;
legend('x[n]','h[n]');
xlabel('n');

%% Assignment 2: Fade-in and -out of convolution
%a
%the length should be equal to N + M -1

%b sample n =[1:8] have no fade-in and fade-out



%% Assignment 4: Frequency response
figure;
nh = [-1 0 1];
h = [ 0.5 1 0.5];
y = fft(h);
m = abs(y);                               % Magnitude
p = unwrap(angle(y));                     % Phase

f = (0:length(y)-1)*1/length(y);        % Frequency vector

subplot(2,1,1)
plot(f,m)
title('Magnitude')
ax = gca;
ax.XTick = [15 40 60 85];

subplot(2,1,2)
plot(f,p*180/pi)
title('Phase')
ax = gca;
ax.XTick = [15 40 60 85];