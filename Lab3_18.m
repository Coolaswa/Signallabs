clear all;
close all;

%18a. 

%b. 

for i = 1:5
    for n = 1:1000
        x(n) = normrnd(0,1);
    end
    u(i) = mean(x);
end