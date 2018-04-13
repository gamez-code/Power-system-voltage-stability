function [Pg1, Qg1, DS,Pperd_l] = PowerBalance (R, G, B, theta1, theta2, V1, V2,Pref)


%---------------------------- Ecuaciones ----------------------------------
%                     ------  Ecuacion 1: ------
%cosenos y senos:
c12 = cos( theta1 - theta2 );
s12 = sin( theta1 - theta2 );

G11V1 = G(1,1) * V1^2;

P12 =  V1 * V2 * ( G(1,2) * c12 + B(1,2) * s12 );
SP1k = P12;

Pn1 = G11V1 + SP1k;
Pg1 = Pn1;

G12V1= (-G(1,2)) * V1^2;

V1V2GB= V1*V2* (G(1,2)*c12 + B(1,2)* s12);


P_12 = G12V1 +  V1V2GB;

DP1 = Pg1 - P_12;

%                     ------  Ecuacion 2: ------
%cosenos y senos:
c21 = cos( theta2 - theta1 );
s21 = sin( theta2 - theta1 );

G21V2= (-G(2,1)) * V2^2;

V2V1GB= V2*V1* (G(2,1)*c21 + B(2,1)* s21);


P_21 = G21V2 +  V2V1GB;

DP2 = -Pref - P_21;

%                   ------ Ecuacion 3: ------
B11V1 = B(1,1) * V1^2;

Q12 =  V1 * V2 * ( G(1,2) * s12 - B(1,2) * c12 );
SQ1k =  Q12;

Qn1 = -B11V1 + SQ1k;
Qg1 = Qn1;


b12o=R(1,5)/2;

B12V1= (B(1,2)-b12o) * V1^2;

V1V2BG= V1*V2* (G(1,2)*s12 - B(1,2)* c12);


Q_12 = B12V1 +  V1V2BG;


DQ1 = Qg1 - Q_12;

%                     ------  Ecuacion 4: ------
%cosenos y senos:
B21V2= (B(2,1)-b12o) * V2^2;

V2V1BG= V2*V1* (G(2,1)*s21 - B(2,1)* c21);


Q_21 = B21V2 +  V2V1BG;

DQ2 = - Q_21;

DS=[DP1+1j*DQ1; DP2+1j*DQ2];
Pperd_l=P_12-P_21;