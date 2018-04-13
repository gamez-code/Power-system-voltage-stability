function [SE] = Ecuacion (x,  y, h, ws, H1,ka, Pm1, ME1, xg1, Iiny, Ybf,...
    Pe1o, Pt_o, Pref,tau,Ybn,Ybo)
%--------------------------- Incognitas -----------------------------------
    delt1 =  x(1);
    omeg1 =  x(2);
    Gt    =  x(3);
    
    delt1o =  y(1);
    omeg1o =  y(2);
    Gto    =  y(3);
    
    [Ybff]=BuiltYbf(Ybf,Ybn,Ybo,Gt);
    [Pe1,Pt]= Algoritm (x,ME1, xg1, Iiny, Ybff);

%--------------------------- Ecuaciones -----------------------------------

    F1 =  - omeg1 + omeg1o + h / ( 4 * H1 ) * ...
        ( Pm1 - Pe1 - ka * omeg1 + Pm1 - Pe1o - ka * omeg1o );
     
    F2 =  - delt1 + delt1o + h * ws/ 2 *( omeg1 + omeg1o );

    F3 =  - Gt + Gto + h / ( 2  * tau ) * ( Pref - Pt_o + Pref  - Pt );
    
    % Sistemas de ecuaciones
    SE =  [ F1, F2, F3];
