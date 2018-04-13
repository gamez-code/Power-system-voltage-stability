function [F,Sd3] = PowerFlow (x, G, B, theta1, V1,Pref)
%---------------------------- Incognitas ----------------------------------

theta2 =  x(1);
V2 =  x(2);
Gt= x(3);

%---------------------------- Ecuaciones ----------------------------------
%                     ------  Ecuacion 1: ------
%cosenos y senos:
c21 = cos( theta2 - theta1 );
s21 = sin( theta2 - theta1 );

G22V2 = G(2,2) * V2^2;

P21 =  V2 * V1 * ( G(2,1) * c21 + B(2,1) * s21 );
SP2k = P21;

Pn2 = G22V2 + SP2k;
Pd2= Pref;
DP2 = -Pd2 - Pn2;

f1 = DP2;
%                   ------ Ecuacion 2: ------
B22V2= B(2,2)*V2^2;

Q21= V1*V2*(G(1,2)*s21-B(1,2)*c21);
SQ2k= Q21;

Qn2=-B22V2+SQ2k;
DQ2= -Qn2;
f2=DQ2;
%                   ------ Ecuacion 3: ------
%cosenos y senos:
f3= Pref-V2^2*Gt;
%----------------- Sistemas de ecuaciones a resolver ----------------------
F=[f1 f2 f3];

