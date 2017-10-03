clear all; close all;

%% Assignment 11: The DFT of a finite length discrete-time signal
figure('NumberTitle', 'off', 'Name', 'Assignment 11: The DFT of a finite length discrete-time signal');
theta = 0:0.001:2*pi;
A = 1; K1 = 11/2; K2 = 1/2; K3 = 0;
X = A.*sin(K1.*theta)./sin(K2.*theta).*exp(-1j*K3.*theta);
plot(theta,abs(X));
grid on;
hold on;
N = 11;
xn = ones(1,11);
Xn = fft(xn,N);
n = 0:2*pi/N:(N-1)/N*2*pi;
stem(n,abs(Xn));

%% Assignment 12: Spectrum of a sine wave of different frequencies
clear all;
figure('NumberTitle', 'off', 'Name', 'Assignment 12: Spectrum of a sine wave of different frequencies');
N = 32;
f1 = 8; f2 = 9; fs = 64;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
n = 0:(N-1);
x1n = sin(theta1.*n); x2n = sin(theta2.*n);
X1n = fft(x1n); X2n = fft(x2n);
theta = 0:2*pi/N:(N-1)/N*2*pi;
subplot(2,1,1);
stem(theta,X1n);
grid on;
subplot(2,1,2);
stem(theta,X2n);

%This can be prevented by making the sample frequency in such a way that
%all the frequencies in the spectrum are a multiple of the frequency

%% Assignment 13: Approximation of the FTD
figure('NumberTitle', 'off', 'Name', 'Assignment 13: Approximation of the FTD');
theta = 0:0.001:2*pi;
A = 1; K1 = 11/2; K2 = 1/2; K3 = 0;
X = A.*sin(K1.*theta)./sin(K2.*theta).*exp(-1j*K3.*theta);
plot(theta,abs(X));
grid on;
hold on;
N = 32;
xn = ones(1,11);
Xn = fft(xn,N);
n = 0:2*pi/N:(N-1)/N*2*pi;
stem(n,abs(Xn));