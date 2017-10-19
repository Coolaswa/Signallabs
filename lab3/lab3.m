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
yn = (x1n + (N-1) .* x2n)./N;
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
    
    