close all; clear all;

%% Assignment 1: Convolution
nx = [-3 -2 -1 0 1 2 3];
x = [1 0 1 0 1 0 1]; %start at n=-3, ends at n=3
h = [-0.5 1 -0.5 ];
y = conv(x,h);
ny = linspace(nx(1) - 1, length(y)+nx(1)-1-1,length(y));

figure('NumberTitle', 'off', 'Name', 'Assignment 1: Convolution');
% subplot(2,1,1);
% stem(nx,x);
% hold on;
% xlim([-4 4]);
% ylim([-1.5 1.5]);
% grid on;
% legend('x[n]','h[n]');
% xlabel('n');
% subplot(2,1,2);
stem(ny,y);
xlim([min([nx,ny])-1 max([nx,ny])+1]);
ylim([min([x,y])-0.5 max([x,y])+0.5]);
grid on;
ylabel('y[n] = h[n]*x[n]');
xlabel('n');

%% Assignment 2: Fade-in and -out of convolution
%a
%the length should be equal to N + M -1

%b sample n =[1:8] have no fade-in and fade-out

%% Assignment 3: Causality
%a) L = 5;
clear all;
x = [1 0 1 0 1 0 1];
nx = -3:3;

h = [1 1 1];
nh_causal = 0:2;
nh = -5:-3;
figure('NumberTitle', 'off', 'Name', 'Assignment 3: Causality');
[ny,y] = convcool(nx, x, nh, h);
[ny_causal,y_causal] = convcool(nx, x, nh_causal, h);
subplot(1,2,1);
stem(ny_causal,y_causal);
xlim([min(ny_causal)-1 max(ny_causal)+1]);
ylim([min(y_causal)-0.5 max(y_causal)+0.5]);
grid on;
ylabel('y[n] = h_{causal}[n]*x[n]');
xlabel('n');
subplot(1,2,2);
stem(ny,y);
xlim([min(ny)-1 max(ny)+1]);
ylim([min(y)-0.5 max(y)+0.5]);
grid on;
ylabel('y[n] = h[n]*x[n]');
xlabel('n');

%% Assignment 4: Frequency response
figure('NumberTitle', 'off', 'Name', 'Assignment 4: Frequency Response');
N = 1000;
nh = [-1 0 1];
h = [0.5 1 0.5];
H = FTD(nh,h,N);
plotMagPhase(H,N);

%% Assignment 5: Frequency response of causal filter
% a&b)
figure('NumberTitle', 'off', 'Name', 'Assignment 5: Frequency response of causal filter');
N = 1000;
nhcausal = [0 1 2];
hcausal = [0.5 1 0.5];
Hcausal = FTD(nhcausal, hcausal, N);
plotMagPhase(Hcausal,N);
%c) the only difference is the phase difference

%% Assignment 8: Sampling a sinusoidal signal
figure('NumberTitle', 'off', 'Name', 'Assignment 8: Sampling a sinusoidal signal');
f=3200; A = 1; phi = 0; fs = 4000;
N = 1000*1/f*fs; %%100 periods
n = 0:(N-1);
theta = 2*pi*(f/fs);
x = A*sin(theta.*n+phi);
stem(1/fs.*n,x);
hold on;
t = 0:0.000001:10;
plot(t,A*sin(2*pi*f.*t+phi));
xlim([0 1/f*2]);
grid on;
soundsc(x,fs);

%% Assignment 9: Visualization of 'aliasing' via time domain
clear all;
fs = 0.5;
f1 = 0.2; f2 = 0.7; f3 = 1.2;
phi = pi/4;
t = 0:0.001:1/f1;
n = 0:length(0:1/fs:1/f1)-1;
figure('NumberTitle', 'off', 'Name', 'Assignment 9: Visualization of aliasing via time domain');
plot(t,cos(2*pi.*f1.*t + phi));
hold on;
plot(t,cos(2*pi.*f2.*t + phi));
plot(t,cos(2*pi.*f3.*t + phi));
stem(n*1/fs,cos(2*pi*f1/fs.*n+phi));

