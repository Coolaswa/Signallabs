clear all;
close all;

%a See slide 93

%b 
x = [1,1,1,0.5];
nx = 0:3;

h = [1,0.75,0.5,0.25];
nh = 0:3;

yp = cconv(x,h,4);
 

%c
[ny,y] = convcool(nx,x,nh,h);