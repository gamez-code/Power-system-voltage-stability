clc
clear all
close all
R=0.1;
f=2;
T=1/f;
dt=T*1e-1;
t=0:dt:17*T;
w=2*pi*f;
V=sin(t).*exp(-0);
P=V.^2/R;
figure(1)
plot(P,V)
grid
figure(2)
plot(t,V)
grid