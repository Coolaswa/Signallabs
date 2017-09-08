function [ ny, y ] = convcool(nx, x, nh, h  )
%CONVCOOL -> 
%   Detailed explanation goes here
    y = conv(x,h);
    ny = linspace(nx(1) + nh(1), length(y)+nx(1) + nh(1) - 1,length(y));
%     stem(nx,x);
%     hold on;
%     stem(ny,y);
%     xlim([min([nx,ny])-1 max([nx,ny])+1]);
%     ylim([min([x,y])-0.5 max([x,y])+0.5]);
%     grid on;
%     legend('x[n]','y[n]');
%     xlabel('n');

end

