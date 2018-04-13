function [sol,Pe1o,Pt_o, V1, V2, Ig1, E1, V]=DinamicSolver(t,y,x, h, ws, H1,ka,...
            Pm1, ME1, xg1, Iiny, Yb, Pref,tau,options,Ybn,Ybo)
for k = 1 : numel(t)
    
    x0 = y;
    [Ybff]=BuiltYbf(Yb,Ybn,Ybo,y(3));
    [Pe1o(k),Pt_o(k), Qe1, V1(k), V2(k),Ig1(k),E1(k),V(:,k)]= Algoritm (x,ME1, xg1, Iiny, Ybff);
    
    
    if k == 1
        x = x;
    else
        [x] = fsolve( @(x)Ecuacion( x,  y, h, ws, H1,ka,...
            Pm1, ME1, xg1, Iiny, Yb,Pe1o(k), Pt_o(k), Pref,tau,Ybn,Ybo), x0, options);
    end
    % Solucion
    sol(:,k)= x';
    y = x;
end