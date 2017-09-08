close all; clear all;

%% Assignment 1: Convolution
nx = [-3 -2 -1 0 1 2 3];
x = [1 0 1 0 1 0 1]; %start at n=-3, ends at n=3
h = [-0.5 1 -0.5 ];
y = conv(x,h);
ny = linspace(nx(1) - 1, length(y)+nx(1)-1-1,length(y));
stem(nx,x);
hold on;
stem(ny,y);

%% Assignment 2: fFade-in and -out of convolution
nh