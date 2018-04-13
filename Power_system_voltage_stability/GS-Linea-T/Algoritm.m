function [Pe1,Pt, Qe1, V1, V2, Ig1,E1,V_nodal]= Algoritm (x,ME1, xg1, Iiny, Y)

%--------------------------- Incognitas -----------------------------------
delt1 =  x(1);
omeg1 =  x(2);
Gt   =  x(3);

%%%%%%%%%%%%%%%%%%% ALGORITMO DE CALCULO DE POTENCIA %%%%%%%%%%%%%%%%%%%%%%

E1 =  ME1 * exp ( 1j * delt1 );

I_iny1 =  E1 / ( 1j * xg1 );

Iiny (1) =   I_iny1;
         
V_nodal =  inv( Y ) * Iiny;
V1 =   V_nodal(1);
V2 =   V_nodal(2);


Vxg1 =  E1 - V1;

Ig1 =  Vxg1 / ( 1j * xg1 );

Pe1 =  real( E1 * conj( Ig1 ) );

Qe1 =  imag( E1 * conj( Ig1 ) );

Pt= abs(V2^2)*Gt;




