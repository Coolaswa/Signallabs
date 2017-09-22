close all;
clear all;

% a. Manually calculating results in (1/(pi*n)) * sin((pi/3)*n)

% b. 
n = -5:5;
h_tilde = (1./(pi.*n)) .* sin((pi/3).*n);
h_tilde(6) = 1/(pi);
figure('NumberTitle', 'off', 'Name', 'Assignment 6: Impulse response of a low pass filter');

stem(n,h_tilde);
xlim([-6 6]);
xlabel('n');
ylabel('h˜');
title('Low pass filter modeled impulse response');


%c. 
n = 0:10000;
Fs = 4000;
x = sin((pi/6) .* n) + sin(((5*pi)/6) .* n);
[n,y] = convcool(-5:5, h_tilde, n, x);
soundsc(y, Fs);
%stem(x);