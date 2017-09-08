close all;
%a) L should equal 5, so that the h[n] equals 1 when n=0,1,2 and equals
%zero elsewhere.

%b) 


x = [1 0 1 0 1 0 1];
nx = -3:3;

h = [1 1 1];
nh_causal = 0:2;
nh = -5:-3;

[ny,y] = convcool(nx, x, nh, h);
[ny_causal,y_causal] = convcool(nx, x, nh_causal, h);

plotcool(ny, y, ny_causal, y_causal);