function [dfidpc,dfidqc,dfidv2]=grad(R, G, B, V1, V2)
theta1=angle(V1);
theta2=angle(V2);
V2=abs(V2);



%---------------------------- Ecuaciones ----------------------------------
%                     ------  Ecuacion 1: ------
%cosenos y senos:
c21 = cos( theta2 - theta1 );
s21 = sin( theta2 - theta1 );

G21V2= 2*(-G(2,1)) * V2;

V2V1GB= (G(2,1)*c21 + B(2,1)* s21);



%                     ------  Ecuacion 2: ------
%cosenos y senos:

b12o=R(1,5)/2;

B21V2= 2*(B(2,1)-b12o) * V2;

V2V1BG= (G(2,1)*s21 - B(2,1)* c21);



%                     ------  Salidas: ------

dfidpc=-V2V1BG;
dfidqc=V2V1GB;
dfidv2=G21V2.*V2V1BG+B21V2.*V2V1GB;
