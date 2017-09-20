function [ output_args ] = plotMagPhase( freqResponse,N )
%PLOTMAGPHASE Summary of this function goes here
%   Detailed explanation goes here
    theta = -pi:2*pi/(N-1):pi;
    m = abs(freqResponse);                               % Magnitude
    p = unwrap(angle(freqResponse));                     % Phase
    figure(9);
    subplot(2,1,1);
    hold on;
    plot(theta,m)
    title('Magnitude');
    ylabel('Magnitude [-]');
    ax = gca;
    ax.XTick = [-pi,-pi/2,0,pi/2,pi];
    grid on;

    subplot(2,1,2)
    plot(theta,p)
    title('Phase')
    ylabel('Phase [rad]');
    ax = gca;
    ax.XTick = [-pi,-pi/2,0,pi/2,pi];
    grid on;

end

