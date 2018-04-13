function [Pg1, Qg1, Qg2, Qg3, DS] = PowerBalance (R, G, B, theta1, theta2, ...
    theta3, V1, V2, V3, Pg2, Sd3)


%---------------------------- Ecuaciones ----------------------------------
%                     ------  Ecuacion 1: ------
%cosenos y senos:
c12 = cos( theta1 - theta2 );
s12 = sin( theta1 - theta2 );
c13 = cos( theta1 - theta3 );
s13 = sin( theta1 - theta3 );

G11V1 = G(1,1) * V1^2;

P12 =  V1 * V2 * ( G(1,2) * c12 + B(1,2) * s12 );
P13 =  V1 * V3 * ( G(1,3) * c13 + B(1,3) * s13 );
SP1k = P13 + P12;

Pn1 = G11V1 + SP1k;
Pg1 = Pn1;

G12V1= (-G(1,2)) * V1^2;
G13V1= (-G(1,3)) * V1^2;

V1V2GB= V1*V2* (G(1,2)*c12 + B(1,2)* s12);
V1V3GB= V1*V3* (G(1,3)*c13 + B(1,3)* s13);


P_12 = G12V1 +  V1V2GB;
P_13 = G13V1 +  V1V3GB;

DP1 = Pg1 - P_12 - P_13;

%                     ------  Ecuacion 2: ------
%cosenos y senos:
c23 = cos( theta2 - theta3 );
s23 = sin( theta2 - theta3 );
c21 = cos( theta2 - theta1 );
s21 = sin( theta2 - theta1 );

G21V2= (-G(2,1)) * V2^2;
G23V2= (-G(2,3)) * V2^2;

V2V1GB= V2*V1* (G(2,1)*c21 + B(2,1)* s21);
V2V3GB= V2*V3* (G(2,3)*c23 + B(2,3)* s23);


P_21 = G21V2 +  V2V1GB;
P_23 = G23V2 +  V2V3GB;

DP2 = Pg2 - P_21 - P_23;

%                   ------ Ecuacion 3: ------
%cosenos y senos:
c32= cos(theta3-theta2);
s32= sin(theta3-theta2);
c31= cos(theta3-theta1);
s31= sin(theta3-theta1);


G31V3= (-G(3,1)) * V3^2;
G32V3= (-G(3,2)) * V3^2;

V3V1GB= V3*V1* (G(3,1)*c31 + B(3,1)* s31);
V3V2GB= V3*V2* (G(3,2)*c32 + B(3,2)* s32);


P_31 = G31V3 +  V3V1GB;
P_32 = G32V3 +  V3V2GB;

Pd3= real(Sd3);
DP3 = -Pd3 - P_31 - P_32;

%                     ------  Ecuacion 4: ------
%cosenos y senos:
B11V1 = B(1,1) * V1^2;

Q12 =  V1 * V2 * ( G(1,2) * s12 - B(1,2) * c12 );
Q13 =  V1 * V3 * ( G(1,3) * s13 - B(1,3) * c13 );
SQ1k = Q13 + Q12;

Qn1 = -B11V1 + SQ1k;
Qg1 = Qn1;


b12o=R(1,5)/2;
b13o=R(2,5)/2;
b23o=R(3,5)/2;

B12V1= (B(1,2)-b12o) * V1^2;
B13V1= (B(1,3)-b13o) * V1^2;

V1V2BG= V1*V2* (G(1,2)*s12 - B(1,2)* c12);
V1V3BG= V1*V3* (G(1,3)*s13 - B(1,3)* c13);


Q_12 = B12V1 +  V1V2BG;

Q_13 = B13V1 +  V1V3BG;

DQ1 = Qg1 - Q_12 - Q_13;

%                     ------  Ecuacion 5: ------
%cosenos y senos:
B22V2 = B(2,2) * V2^2;

Q21 =  V2 * V1 * ( G(2,1) * s21 - B(2,1) * c21 );
Q23 =  V2 * V3 * ( G(2,3) * s23 - B(2,3) * c23 );
SQ2k = Q23 + Q21;

Qn2 = -B22V2 + SQ2k;
Qg2 = Qn2;


B21V2= (B(2,1)-b12o) * V2^2;
B23V2= (B(2,3)-b23o) * V2^2;

V2V1BG= V2*V1* (G(2,1)*s21 - B(2,1)* c21);
V2V3BG= V2*V3* (G(2,3)*s23 - B(2,3)* c23);


Q_21 = B21V2 +  V2V1BG;
Q_23 = B23V2 +  V2V3BG;

DQ2 = Qg2 - Q_21 - Q_23;

%                   ------ Ecuacion 6: ------
%cosenos y senos:



B31V3= (B(3,1)-b13o) * V3^2;
B32V3= (B(3,2)-b23o) * V3^2;

V3V1BG= V3*V1* (G(3,1)*s31 - B(3,1)* c31);
V3V2BG= V3*V2* (G(3,2)*s32 - B(3,2)* c32);


Q_31 = B31V3 +  V3V1BG;
Q_32 = B32V3 +  V3V2BG;

wC=inv(-1j*1.5);
Qd3= imag(Sd3);
Qg3 = V3^2*abs(wC);
DQ3 = Qg3 - Qd3 - Q_31 - Q_32;

DS=[DP1+1j*DQ1; DP2+1j*DQ2; DP3+1j*DQ3 ];