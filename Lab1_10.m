clear all; close all;

[y, fsin] = audioread('scale.wav');
K=2;
N =51;
nh = -(N-1)/2:(N-1)/2;
theta_c = pi/2;
h = theta_c/pi*sinc(nh*theta_c/pi);

%% downsampling
%a)
audiowrite('out_dec1.wav',y./max(y),fsin);

%b)
y_out = y(1:2:length(y));
audiowrite('out_dec2.wav',y_out./max(y_out),fsin/2);

%c) 
[~, y_out] = convcool(0:length(y),y,nh,h);
y_out = y_out(1:2:length(y_out));
audiowrite('out_dec3.wav',y_out./max(y_out),fsin/2);

%% upsampling
%a)
audiowrite('out_up1.wav',y./max(y),fsin);

%b)
y_out(1:2:(2*length(y))) = y;
audiowrite('out_up2.wav',y_out./max(y_out),fsin*2);

%c)
y_out(1:2:(2*length(y))) = y;
[~, y_out] = convcool(0:length(y_out),y_out,nh,h);
audiowrite('out_up3.wav',y_out./max(y_out),fsin*2);
