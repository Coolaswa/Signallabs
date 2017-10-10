clear all; close all;

%% Assignment 11: The DFT of a finite length discrete-time signal
figure('NumberTitle', 'off', 'Name', 'Assignment 11: The DFT of a finite length discrete-time signal');
theta = -pi:0.001:pi;
A = 1; K1 = 11/2; K2 = 1/2; K3 = 0;
X = A.*sin(K1.*theta)./sin(K2.*theta).*exp(-1j*K3.*theta);
plot(theta,abs(X));
grid on;
hold on;
Na = 11;
xn = ones(1,11);
Xn = fft(xn,Na);
n = (-(Na-1)/2:(Na-1)/2)*2*pi/Na;
stem(n,fftshift(abs(Xn)));
legend('abs(X(e^{j0}))', 'abs(X[k])');
xlabel('\theta [rad/s]');
ylabel('|X(\theta|');

%% Assignment 12: Spectrum of a sine wave of different frequencies
clear all;
figure('NumberTitle', 'off', 'Name', 'Assignment 12: Spectrum of a sine wave of different frequencies');
Na = 32;
f1 = 8; f2 = 9; fs = 64;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
n = 0:(Na-1);
x1n = sin(theta1.*n); x2n = sin(theta2.*n);
X1n = fft(x1n); X2n = fft(x2n);
theta = -pi:2*pi/Na:(Na-1)/Na*pi;
subplot(2,1,1);
stem(theta,abs(X1n));
title('DFT of x_{1}[n]');
grid on;
xlim([-3.5,3.5]);
xlabel('\theta [rad/s]');
ylabel('X_1(\theta)');
subplot(2,1,2);
stem(theta,abs(X2n));
grid on;
title('DFT of x_{2}[n]');
xlim([-3.5,3.5]);
xlabel('\theta [rad/s]');
ylabel('X_1(\theta)');

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
xlabel('\theta [rad/s]');
ylabel('|X(\theta)|');
legend('abs(X(e^{j0}))', 'abs(X[k])');
title('Approximation of the FTD');

%% Assignment 14: Calculating the minimum resolution of a spectrumclear all;
%b) 
figure('NumberTitle', 'off', 'Name', 'Assignment 14: Calculating the minimum resolution of a spectrum');
f1 = 175; f2 = 200; fs = 1000;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
Na = ceil(0.89*2*pi/abs(theta1-theta2));
n = 0:Na-1;
xn = sin(theta1.*n) + sin(theta2.*n);
Nplot = 5*Na;
Xn = fft(xn,Nplot);
theta = -pi:2*pi/Nplot:(pi-1/Nplot*2*pi);%0:2*pi/Nplot:((Nplot)-1)/(Nplot)*2*pi;
%theta(length(theta)/2:end) = theta(length(theta)/2:end) - 2*pi;
subplot(2,3,1);
plot(theta,fftshift(abs(Xn)));
title('b');
xlabel('\theta [rad/s]');
ylabel('|W(\theta) \ast X(\theta)|');
grid on;
%c)
f1 = 185; f2 = 200; fs = 1000;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
subplot(2,3,2);
xn = sin(theta1.*n) + sin(theta2.*n);
Nplot = 10*Na;
Xn = fft(xn,Nplot);
theta = -pi:2*pi/Nplot:(pi-1/Nplot*2*pi);%0:2*pi/Nplot:((5*Nplot)-1)/(5*Nplot)*2*pi;
plot(theta,fftshift(abs(Xn)));
title('c')
xlabel('\theta [rad/s]');
ylabel('|W(\theta) \ast X(\theta)|');
grid on;
%d)
f1 = 185; f2 = 200; fs = 1000;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
Nb = ceil(0.89*2*pi/abs(theta1-theta2));
n = 0:Nb-1;
xn = sin(theta1.*n) + sin(theta2.*n);
Nplot = 1000;
Xn = fft(xn,Nplot);
theta = -pi:2*pi/Nplot:(pi-1/Nplot*2*pi);
subplot(2,3,3);
plot(theta,fftshift(abs(Xn)));
title('d');
xlabel('\theta [rad/s]');
ylabel('|W(\theta) \ast X(\theta)|');
grid on;

%e)
f1 = 185; f2 = 200; fs = 1000;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
Nc = ceil(1.44*2*pi/abs(theta1-theta2));
n = 0:Nc-1;
wn = (1-cos(2*pi/Nc.*n))/2;
xn = sin(theta1.*n) + sin(theta2.*n);
xn = wn.*xn;
Nplot = 5*Nc;
Xn = fft(xn,Nplot);
theta = -pi:2*pi/Nplot:(pi-1/Nplot*2*pi);%theta = 0:2*pi/Nplot:((5*Nplot)-1)/(5*Nplot)*2*pi;
subplot(2,3,4);
plot(theta,fftshift(abs(Xn)));
title('e');
xlabel('\theta [rad/s]');
ylabel('|W(\theta) \ast X(\theta)|');
grid on;

%f)
f1 = 185; f2 = 200; fs = 1000;
theta1 = 2*pi*f1/fs; theta2 = 2*pi*f2/fs;
Nc = ceil(1.44*2*pi/abs(theta1-theta2));
n = 0:Nc-1;
wn = (1-cos(2*pi/Nc.*n));
xn = sin(theta1.*n) + 0.35*sin(theta2.*n);
xn = wn.*xn;
Nplot = 5*Nc;
Xn = fft(xn,Nplot);
theta = -pi:2*pi/Nplot:(pi-1/Nplot*2*pi);%theta = 0:2*pi/Nplot:((5*Nplot)-1)/(5*Nplot)*2*pi;
subplot(2,3,5);
plot(theta,fftshift(abs(Xn)));
xlabel('\theta [rad/s]');
ylabel('|W(\theta) \ast X(\theta)|');
grid on;
title('e');

%% Assignment 16: Frequency plots of LPF
clear all;

for N = [5 11 101]
    n = -(N-1)/2:1:(N-1)/2;
    hn = 0.5 * sinc(pi.*n/2);
    figure('NumberTitle', 'off', 'Name', sprintf('Assignment 16, N = %d', N));
    freqz(hn);
    title(sprintf('Spectrum of LPF of size %d',N));
end

%% Assignment 17: Filter design for audio signals 1
clear all;
[y, fs] = audioread('audio_sin.wav ');
Y = fft(y);
N = length(y);
thetan = -pi:2*pi/(N-1):pi;
figure('NumberTitle', 'off', 'Name', 'Assignment 17');
plot(thetan/2/pi*fs,abs(fftshift(Y))); %fa = 8820Hz
xlabel('frequency [Hz]');
ylabel('|X(\theta)|');
grid on;
load('bandstopfilter.mat');
yn = conv(y,Num);
audiowrite('new.wav',yn./max(yn),fs);