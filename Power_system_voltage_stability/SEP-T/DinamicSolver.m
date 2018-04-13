function [sol,Pe1o, Pe2o, Pe3o, Qe1o, Qe2o, V_nodal, ...
         I3o,I2o,I1o]=DinamicSolver(t,y,x, h, ws, H1,H2,ka,...
            Pm1,Pm2, ME1, ME2, xg1, xg2, Iiny, Yb, Pref,tau,options,R,G,B)
for k = 1 : numel(t)

    x0 = y;
    [Pe1o(k), Pe2o(k), Pe3o(k), Qe1o(k), Qe2o(k), V_nodal(:,k), ...
        Q1, Q2, I3o(k),I2o(k),I1o(k) ] = Algoritm (y, ME1, ME2, xg1, xg2, Iiny, Yb );
       
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( @(x)Ecuacion(x,  y, h, ws, H1, ...
            H2, ka, Pm1, Pm2,tau, Pref, ME1, ME2, xg1, xg2,...
            Iiny, Yb, Pe1o(k), Pe2o(k), Pe3o(k)), x0, options);
    end
    % Solucion
    sol(:,k) = x';
    y = x;
    
    if numel(Iiny)==3
    [DS] = PowerBalance2 (R, G, B, Pe1o(k), Pe2o(k), Pe3o(k), Q1,...
        Q2, V_nodal(:,k));
    if DS > 1e-10
        error('falla en el balance de potencia')
    end
    elseif numel(Iiny)==4
         [DS] = PowerBalance3 (R, G, B, Pe1o(k), Pe2o(k), Pe3o(k), Q1, Q2,...
             V_nodal(:,k));
    if sum(abs(DS)) > 1e-10
       error('falla en el balance de potencia')
    end
    end
end
end