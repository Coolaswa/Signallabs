N = 1000;
samp = 10;
nh = 0:samp;
theta_c = 0.63/4;
h = 1/samp*ones(1,samp);
h = 2*theta_c/pi*sinc(nh.*theta_c/pi);
H = FTD(nh,h,N);
plotMagPhase(H,N);