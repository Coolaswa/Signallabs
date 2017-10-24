clear all;
close all;

%18a. 

%b. 

for n = 1:1000
    for i = 1:5
        x(i,n) = normrnd(0,1);
    end
    u1(n) = x(1,n);
    u2(n) = mean(x(1:2,n));
    u5(n) = mean(x(:,n));
end
fig = figure;
subplot(3,1,1);
plot(u1);
ylim([-4 4]);
title('N = 1: 1000 random samples');
xlabel('Samples');
ylabel('Amplitude');
subplot(3,1,2);
plot(u2);
ylim([-4 4]);
title('N = 2: 2 random samples averaged');
xlabel('Samples');
ylabel('Amplitude');
% hold on;
% for n = 1:1000
%     test(n) = normrnd(0,sqrt(0.5));
% end
% plot(test);
subplot(3,1,3);
plot(u5);
ylim([-4 4]);
title('N = 5: 5 random samples averaged');
xlabel('Samples');
ylabel('Amplitude');
% hold on;
% for n = 1:1000
%     test(n) = normrnd(0,sqrt(0.2));
% end
% plot(test);
saveas(fig,'Assignment18.png');