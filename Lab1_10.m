clear all; close all;

[y, fsin] = audioread('scale.wav');
K=2;
N = 11;
nh = -(N-1)/2:(N-1)/2;
theta_c = pi/K;
h = theta_c/pi*sinc(nh*theta_c/pi);

%% downsampling
%a)
audiowrite('out_dec1.wav',y,fsin);

%b)
y_out = y(1:2:length(y));
audiowrite('out_dec2.wav',y_out,fsin/2);

%c) 
[~, y_out] = convcool(0:length(y),y,nh,h);
y_out = y_out(1:2:length(y_out));
audiowrite('out_dec3.wav',y_out,fsin/2);
