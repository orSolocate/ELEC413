clear all;
clc;

mat100130=load("orbahari_100130_1082 (1).mat");
mat10090=load("orbahari_10090_1085.mat");
mat100170=load("orbahari_100170_1083 (1).mat");
mat100210=load("orbahari_100210_1084 (1).mat");
%%
figure(82);
plot(mat100130.scandata.wavelength*1e9,mat100130.scandata.power(:,2));
%R=mat100130.scandata.power(:,1);
Rpeak100130=min(mat100130.scandata.power(2200:3000,1))
loss100130=10*log10(1-10^(Rpeak100130/10))
%%
figure(83);
plot(mat100170.scandata.wavelength*1e9,mat100170.scandata.power(:,2));
Rpeak100170=min(mat100170.scandata.power(2200:3000,1))
loss100170=10*log10(1-10^(Rpeak100170/10))
%%
figure(84);
plot(mat100210.scandata.wavelength*1e9,mat100210.scandata.power(:,2));
Rpeak100210=min(mat100210.scandata.power(2200:3000,1))
loss100210=10*log10(1-10^(Rpeak100210/10))
%%
figure(85);
plot(mat10090.scandata.wavelength*1e9,mat10090.scandata.power(:,2));
Rpeak10090=min(mat10090.scandata.power(2200:3000,1))
loss10090=10*log10(1-10^(Rpeak10090/10))

%% loading unused results
mat1087=load("orbahari_60130_1087 (1).mat");
mat1080=load("orbahari_60170_1080 (1).mat");
mat1081=load("orbahari_60210_1081 (1).mat");
mat1086=load("orbahari_140210_1086.mat");
tester=load("orbahari_tester_1088.mat");

%%
plot(tester.scandata.wavelength*1e9,tester.scandata.power(:,2));
title("Tester Results")
xlabel('wavelength [nm]')
ylabel('poewr [dB]')
hold on
plot(tester.scandata.wavelength*1e9,tester.scandata.power(:,1));
legend("T","R")
hold off;