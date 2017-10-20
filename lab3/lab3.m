clear all; close all;

%% Assignment 19: Mathematical expressions for variance and correlation coefficient
% a)
N = 1:100;
sigmax2 = (2-2*N+N.^2)./N.^2;
plot(N,sigmax2);
xlabel('N [-]');
ylabel('\sigma_{yN}^2');

% b)
rho = (N-1)./N./sqrt(sigmax2);
figure;
plot(N,rho);

%% Assignment 20: Scatter plots and empirically evaluate correlation coefficient
% a)
x1n = normrnd(0,1,1000,1);
x2n = normrnd(0,1,1000,1);
figure;
hold on;
for N = [1,2,5,100]
    yn = (x1n + (N-1) .* x2n)./N;
    plot(x2n,yn,'.');
end
legend('N = 1','N = 2', 'N = 5', 'N = 100');

% b)
figure;
rhohat = [];
for N = 1:100
    yn = (x1n + (N-1) .* x2n)./N;
    rhohat = [rhohat (sum((x2n - mean(x2n)).*(yn-mean(yn))))./sqrt(sum((x2n-mean(x2n)).^2).*sum((yn-mean(yn)).^2))];
end;
plot(1:100,rhohat);

% c) The different samples of the scatter plots are more on one line when
% the normalized cross correlation coefficient is higher.

%% Assignment 21: Correlation function and power spectral density (PSD) function
%zie schrift


%% Assignment 22: Scatter plots
n = 1:900;
N = 10;
xn = normrnd(0,1,1000,1);
yn = 1/3.*([xn; 0; 0] + [0;xn;0] + [0;0;xn]);
figure;
for i = 1:4
    subplot(2,2,i);
    plot(yn(n),yn(n+i),'.');
end

%% Assignment 23: Empirical correlation and PSD function
L = 11;
M = length(yn);
ry = [];
for l = -(L-1):1:(L-1)
    rytemp = 0;
    for n = 1:(M-abs(l))
        rytemp = rytemp + 1./M.*yn(n).*yn(n+abs(l));
    end
    ry = [ry rytemp];
end
l = -(L-1):1:(L-1);
figure;
stem(l,ry);
hold on;
rytheory = 1/9.*[zeros(1,8) 1 2 3 2 1 zeros(1,8)];
stem(l,rytheory);

% b) for l = 0, there is a maximum

% c)
Pxest = 0;
l = -(L-1):(L-1);
theta = 0:0.001:2*pi;
for i = 1:length(l)
    Pxest = Pxest + ry(i).*exp(-1j.*theta.*l(i));    
end
figure;
plot(theta,abs(Pxest));
hold on;
plot(theta,1/9.*(3+4.*cos(theta)+2.*cos(2.*theta)));
legend('estimation', 'theory');

%% Assignment 24: Expression for cross-correlation
%a) zie schrift
%b) You can estimate the delay by measuring the signal and then calculating
%the estimation of the autocorrelation.

%% Assignment 25: Test cross-correlation function
clear all;
xn = normrnd(0,1,1000,1);
zn = normrnd(0,1,1000,1);
yn = xn + zn;
l = -10:10;
rxy = crosscor(xn,yn,l);
figure;
stem(l,rxy);

%% Assignment 26: Estimate delay for radar data
%a);
load('radar.mat');
n = 1:length(trans);
plot(n,trans,n,received);

%b)
rx = [];
for l = -100:100;
    rx =[rx crosscor2(trans,trans,l)];
end
rxy = [];
for l = -100:100
    rxy = [rxy crosscor2(trans,received,l)];
end
l = -100:100;
figure;
stem(l,rx);
hold on;
stem(l,rxy);

function ry = crosscor(x,y,l)
    ry = 0;
    M = length(x);
    for n = 1:(M-abs(l))
        ry = ry + 1./M.*x(n).*y(n+abs(l));
    end
end

function ry = crosscor2(x,y,l)
    ry = 0;
    M = length(x);
    if l > 0
        for n = 1:(M-(l))
            ry = ry + 1./M.*x(n).*y(n+(l));
        end
    else
        for n = (abs(l)+1):(M)
            ry = ry + 1/M.*x(n).*y(n+l);
        end
    end
end



    
    