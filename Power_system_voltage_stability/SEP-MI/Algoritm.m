function [Pe1, Pe2, Pe3, Qe1, Qe2, V_nodal, Q1, Q2, Peje, Pm, Qe3, Pf,I3,Veje ]...
    = Algoritm (x,ME1, ME2, xg1, xg2, Iiny, Yb, Ze, Zm, Zr, kc, m, kf)

%--------------------------- Incognitas -----------------------------------
delt1 =  x(1);
delt2 =  x(2);
omeg1 =  x(3);
omeg2 =  x(4);
s    =  x(5);

%-----------------Calculo de la Zeq de la máquina de indución--------------
[Peje, Pm, Pf, Sd3,Pd3,Qd3,In,Zeq,Veje]=MI(Zr, Zm, Ze,0,0, s,0,0,0);
Yeq= inv(Zeq);

%%%%%%%%%%%%%%%%%%% ALGORITMO DE CALCULO DE POTENCIA %%%%%%%%%%%%%%%%%%%%%%

E1 =  ME1 * exp ( 1j * delt1 );
E2 =  ME2 * exp ( 1j * delt2 );

I_iny1 =  E1 / ( 1j * xg1 );
I_iny2 =  E2 / ( 1j * xg2 );

Iiny (1) =   I_iny1;
Iiny (2) =   I_iny2; 

Yb(3,3) = Yb(3,3)+Yeq;
Y = Yb;
         
V_nodal =  inv( Y ) * Iiny;
V1 =   V_nodal(1);
V2 =   V_nodal(2);
V3 =   V_nodal(3);

Vxg1 =  E1 - V1;
Vxg2 =  E2 - V2;

Ig1 =  Vxg1 / ( 1j * xg1 );
Ig2 =  Vxg2 / ( 1j * xg2 );

Pe1 =  real( E1 * conj( Ig1 ) );
Pe2 =  real( E2 * conj( Ig2 ) );
Qe1 =  imag( E1 * conj( Ig1 ) );
Qe2 =  imag( E2 * conj( Ig2 ) );

[Peje, Pm, Pf, Sd3,Pe3,Qe3,I3,Zeq,Veje]=MI(Zr, Zm, Ze,abs(V3),angle(V3), s,kc,m,kf); 

Q1 =  imag( V1 * conj( Ig1 ) );
Q2 =  imag( V2 * conj( Ig2 ) );




