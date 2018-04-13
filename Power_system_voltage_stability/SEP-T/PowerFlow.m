function F = PowerFlow (x, G, B, theta1, V1, V2, Pg2, Pref)
%---------------------------- Incognitas ----------------------------------

theta2 =  x(1);
theta3 =  x(2);
V3 =  x(3);
Gc =  x(4);

%---------------------------- Ecuaciones ----------------------------------
%                     ------  Ecuacion 1: ------
%cosenos y senos:
c23 = cos( theta2 - theta3 );
s23 = sin( theta2 - theta3 );
c21 = cos( theta2 - theta1 );
s21 = sin( theta2 - theta1 );

G22V2 = G(2,2) * V2^2;

P23 =  V2 * V3 * ( G(2,3) * c23 + B(2,3) * s23 );
P21 =  V2 * V1 * ( G(2,1) * c21 + B(2,1) * s21 );
SP2k = P21 + P23;

Pn2 = G22V2 + SP2k;
DP2 = Pg2 - Pn2;

f1 = DP2;
%                   ------ Ecuacion 2: ------
%cosenos y senos:
c32= cos(theta3-theta2);
s32= sin(theta3-theta2);
c31= cos(theta3-theta1);
s31= sin(theta3-theta1);

G33V3= G(3,3)*V3^2;

P32= V3*V2*(G(3,2)*c32+B(3,2)*s32);
P31= V3*V1*(G(3,1)*c31+B(3,1)*s31);
SP3k= P31+P32;

Pn3= G33V3+SP3k;
DP3= -Pref-Pn3;

f2= DP3;

%                   ------ Ecuacion 3: ------
%cosenos y senos:

B33V3= B(3,3)*V3^2;

Q32= V3*V2*(G(3,2)*s32-B(3,2)*c32);
Q31= V3*V1*(G(3,1)*s31-B(3,1)*c31);
SQ3k= Q31+Q32;

Qn3=-B33V3+SQ3k;
DQ3= -Qn3;

f3= DQ3;

%                   ------ Ecuacion 4: ------
% Conductancia:
f4= -Pref+V3^2*Gc;

%----------------- Sistemas de ecuaciones a resolver ----------------------
F=[f1 f2 f3 f4];

