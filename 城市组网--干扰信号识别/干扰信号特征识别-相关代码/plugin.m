%%
clear
clc
close

fs = 3000;
t = 0:1/fs:1-1/fs;

x1 = chirp(t,300,t(end),1300,'quadratic',0,'convex') +randn(size(t))/100;


x2 = exp(2j*pi*100*cos(2*pi*2*t)) + randn(size(t))/100;

[p,f,t] = pspectrum(x2,fs,'spectrogram');

waterfall(f,t,p')
xlabel('Frequency (Hz)')
ylabel('Time (seconds)')
wtf = gca;
wtf.XDir = 'reverse';
view([30 45])


%%
clear
close all
clc

load cap.mat
cap_R = R_sorted;


load not_cap.mat

ncap_R = R_sorted;


