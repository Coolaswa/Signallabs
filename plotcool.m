function [ output_args ] = plotcool( nx, x, nx, y )
%PLOTCOOL Summary of this function goes here
%  Detailed explanation goes here
    stem(nx,x);
    hold on;
    stem(ny,y);
    xlim([min([nx,ny])-1 max([nx,ny])+1]);
    ylim([min([x,y])-0.5 max([x,y])+0.5]);
    grid on;
    legend('x[n]','y[n]');
    xlabel('n');

end

