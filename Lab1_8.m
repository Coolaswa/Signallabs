figure('NumberTitle', 'off', 'Name', 'Assignment 8: Sampling a sinusoidal signal');
f=800; A = 1; phi = 0; fs = 4096;
N = 5000*1/f*fs; %%100 periods
n = 0:(N-1);
theta = 2*pi*(f/fs);
x = A*sin(theta.*n+phi);
stem(1/fs.*n,x);
hold on;
t = 0:0.000001:10;
plot(t,A*sin(2*pi*f.*t+phi));
xlim([0 1/f*2]);
legend('x[n*Ts]', 'x_c(t)');
xlabel('Time [s]');
grid on;
sound(x,fs);