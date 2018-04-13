function [SE] = Ecuacion (x,  y, h, ws, H1, H2,...
            ka, Pm1, Pm2, H3, kc, m, ME1, ME2, xg1, xg2, Iiny, Yb, Pe1o,...
            Pe2o, Peje_o, Pm_o, Ze, Zm, Zr, kf, Pf_o)
%--------------------------- Incognitas -----------------------------------
    delt1 =  x(1);
    delt2 =  x(2);
    omeg1 =  x(3);
    omeg2 =  x(4);
    s    =  x(5);
    
    delt1o =  y(1);
    delt2o =  y(2);
    omeg1o =  y(3);
    omeg2o =  y(4);
    so    =  y(5);
    
    [Pe1, Pe2, Pe3, Qe1, Qe2, V, Q1, Q2, Peje, Pm, Qe3, Pf,I3 ] = Algoritm (x, ME1, ME2, xg1,...
        xg2, Iiny, Yb, Ze, Zm, Zr,  kc, m, kf );
%--------------------------- Ecuaciones -----------------------------------

    F1 =  - omeg1 + omeg1o + h / ( 4 * H1 ) * ...
        ( Pm1 - Pe1 - ka * omeg1 + Pm1 - Pe1o - ka * omeg1o );
    
    F2 =  - omeg2 + omeg2o + h  / ( 4 * H2 ) * ...
        ( Pm2 - Pe2 - ka * omeg2 + Pm2 - Pe2o - ka * omeg2o );
    
    F3 =  - delt1 + delt1o + h * ws/ 2 *( omeg1 + omeg1o );
    F4 =  - delt2 + delt2o + h * ws/ 2 *( omeg2 + omeg2o );
    
    F5 =  - s + so + h / ( 4  * H3 ) * ( Pm_o +Pf_o- Peje_o + Pm+Pf - Peje );
    
    % Sistemas de ecuaciones
    SE =  [ F1, F2, F3, F4, F5];
    