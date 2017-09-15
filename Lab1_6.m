close all;
clear all;

% a. Manually calculating results in (1/(2*pi*n)) * sin((pi/3)*n)

% b. 
n = -5:5;
h_tilde = (1./(2*pi.*n)) .* sin((pi/3).*n);
h_tilde(6) = 1/(2*pi);
%stem(n,h_tilde);

%c. 
n = 0:10000;
Fs = 4000;
x = sin((pi/6) .* n) + sin(((5*pi)/6) .* n);
[n,y] = convcool(-5:5, h_tilde, n, x);
soundsc(y, Fs);
%stem(x);