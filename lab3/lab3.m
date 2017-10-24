clear all; close all;

%% Assignment 19: Mathematical expressions for variance and correlation coefficient
% a)
fig = figure;
N = 1:100;
sigmax2 = (2-2*N+N.^2)./N.^2;
plot(N,sigmax2);
xlabel('N [-]');
ylabel('\sigma_{yN}^2');
grid on;
saveas(fig,'Assignment19a.png');

% b)
rho = (N-1)./N./sqrt(sigmax2);
fig = figure;
plot(N,rho);
xlabel('N [-]');
ylabel('\rho_{x2,yN}[0]');
grid on
saveas(fig,'Assignment19b.png');

%% Assignment 20: Scatter plots and empirically evaluate correlation coefficient
% a)
x1n = normrnd(0,1,1000,1);
x2n = normrnd(0,1,1000,1);
fig  = figure;
N = [1,2,5,100];
for i=1:4
    yn = (x1n + (N(i)-1) .* x2n)./N(i);
    subplot(2,2,i);
    plot(x2n,yn,'.');
    xlabel('x_2[n]');
    ylabel('y_N[n]');
    grid on;
    legendInfo = ['N = ' num2str(N(i))];
    legend(legendInfo);
    
end
grid on;
xlabel('x_2[n]');
ylabel('y_N[n]');
saveas(fig,'Assignment20a.png');

% b)
fig = figure;
rhohat = [];
for N = 1:100
    yn = (x1n + (N-1) .* x2n)./N;
    rhohat = [rhohat (sum((x2n - mean(x2n)).*(yn-mean(yn))))./sqrt(sum((x2n-mean(x2n)).^2).*sum((yn-mean(yn)).^2))];
end;
plot(1:100,rhohat);
grid on;
xlabel('N');
N = 1:100;
hold on;
rho = (N-1)./N./sqrt(sigmax2);
plot(N,rho);

legend({'$\hat{\rho}_{x2,yN}[0]$','$\rho_{x2,yN}[0]$'},'Interpreter','latex');
saveas(fig,'Assignment20b.png');

% c) The different samples of the scatter plots are more on one line when
% the normalized cross correlation coefficient is higher.

%% Assignment 21: Correlation function and power spectral density (PSD) function
%zie schrift


%% Assignment 22: Scatter plots
n = 1:900;
xn = normrnd(0,1,1000,1);
yn = 1/3.*([xn; 0; 0] + [0;xn;0] + [0;0;xn]);
fig = figure;
for i = 1:4
    subplot(2,2,i);
    plot(yn(n),yn(n+i),'.');
    xlabel('y_N[n]');
    yl = ['y_N[n + ' num2str(i) ']'];
    ylabel(yl);
    grid on;
end
saveas(fig,'Assignment22.png');

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
fig = figure;
stem(l,ry);
hold on;
rytheory = 1/9.*[zeros(1,8) 1 2 3 2 1 zeros(1,8)];
stem(l,rytheory);
legend({'$\hat{r}_y[l]$','$r_y[l]$'},'Interpreter','latex');
grid on;
xlabel('l');
ylabel('$\hat{r}_y[l]/r_y[l]');
saveas(fig,'Assignment23a.png');
% b) for l = 0, there is a maximum

% c)
Pxest = 0;
l = -(L-1):(L-1);
theta = 0:0.001:2*pi;
for i = 1:length(l)
    Pxest = Pxest + ry(i).*exp(-1j.*theta.*l(i));    
end
fig=figure;
plot(theta,abs(Pxest));
hold on;
plot(theta,1/9.*(3+4.*cos(theta)+2.*cos(2.*theta)));
legend('estimation', 'theory');
grid on;
xlabel('\theta');
ylabel('P_y(e^{j\theta})');
saveas(fig,'Assignment23c.png');

%% Assignment 24: Expression for cross-correlation
%a) zie schrift
%b) You can estimate the delay by measuring the signal and then calculating
%the estimation of the autocorrelation.

%% Assignment 25: Test cross-correlation function
clear all;
xn = normrnd(0,1,1000,1);
zn = normrnd(0,1,1000,1);
yn = xn + zn;
rxy = [];
for l = -10:10
    rxy = [rxy crosscor2(xn,yn,l)];
end
l = -10:10;
fig = figure;
stem(l,rxy);
grid on;
xlabel('l');
ylabel('r_{xy}[l]');
saveas(fig,'Assignment25a.png');

%% Assignment 26: Estimate delay for radar data
%a);
fig = figure;
load('radar.mat');
n = 1:length(trans);
subplot(2,1,1);
plot(n,trans);
xlabel('n');
ylabel('trans[n]');
grid on;
subplot(2,1,2);
plot(n,received);
xlabel('n');
ylabel('received[n]');
grid on;
saveas(fig,'Assignment26a.png');
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
fig = figure;
subplot(2,1,1);
stem(l,rx);
xlabel('l');
ylabel('r_{trans}[l]');
grid on;
hold on;
subplot(2,1,2);
stem(l,rxy);
xlabel('l');
ylabel('r_{trans,received}[l]');
grid on;
saveas(fig,'Assignment26b.png');
figure;
subplot(2,1,1);
test = abs(fft(trans));
plot(test);
subplot(2,1,2);
ryy = [];
for l = -100:100
    ryy = [ryy crosscor2(received,received,l)];
end
stem(ryy);

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



    
    