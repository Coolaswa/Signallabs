clear all; close all;
load('radar.mat');
load('lowpass.mat');
TRANS = fft(trans);
theta=linspace(-pi,pi,length(TRANS));
fig = figure;
plot(theta,fftshift(abs(TRANS)));
grid on;
xlabel('\theta');
ylabel('FTD magnitude of trans');
saveas(fig,'Assignment26d.png');

%b)
rx = [];
rxy = [];
for l = -100:100
    rxy = [rxy crosscor2(trans,received,l)];
end
rxy2=[];
for l = -100:100
    rxy2 = [rxy2 crosscor2(trans,conv(received,lowpass),l)];
end
l = -100:100;
fig = figure;
subplot(2,1,1);
stem(l,rxy);
xlabel('l');
ylabel('r_{trans,received,normal}[l]');
grid on;
hold on;
subplot(2,1,2);
stem(l,rxy2);
xlabel('l');
ylabel('r_{trans,received,lowpass}[l]');
grid on;

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
