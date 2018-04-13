%%%%%%%%%%%%%% ALGORITMO DINAMICO DE ESTABILIDAD DE VOLTAJE %%%%%%%%%%%%%%%
%Prof. Juan Bermudez
%Br. Manuel Jimenez
clc
clear all
close all
tic

options=optimset('Display','off');

% Use la topologia del sistema de los apuntes JB de
% la Pag 151, con una modificación, la carga no es 
% una potencia constante sino un motor de inducción.
% modelamos la carga con el modelo en régimen permanente
% 
% Corriendo el flujo de carga, Hice una verificación 
% haciendo una balance de potencia.


%%%%%%%%%%%%%%%%%%%%%%%% DATOS DEL PROBLEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Datos del tiempos:-----------------------------------
f=60;                           % frecuencia
ws=2*pi*f;                      % frecuencia angular
T=1/f;                          % Periodo
t_i=100*T;                       % tiempo de inicio de la falla
t_d=30*T;   %tiempo critico 50 ciclos estabilidad de voltaje y angulo                    % tiempo de despeje
%Nota tiempo critico de despeje con solo estabilidad de angulo es de 38
%ciclos, insensible a variaciones del reostato. Constate de tiempo del
%reostato igual a infinito.
t_f=360*T;                      % tiempo final de interacion

h=T/2%T*1e-1;                            % Delta t, pasos de integracion
t1=0:h:t_i;                     % tiempo pre-falla
t2=t_i:h:t_i+t_d;               % tiempo de falla
t3=t_i+t_d:h:t_i+t_d+t_f;       % tiempo post-falla
t=[ t1 t2 t3 ];            
% Vector de tiempo


%--------------------Datos de los generadores------------------------------
%Generador 1:
% Tipo Slack:
V1=1;                           % Voltaje en bornes 
theta1=0;                       % Angulo en bornes
H1=inf%10;                          % Constante de inercia de la mï¿½quina 1
xg1=0.1;                        % Impedancia del generador 1

%Generador 2:
V2=1.01;                        % Voltaje en bornes
H2=inf%5;                           % Constate de inercia se la mï¿½quina 2
Pg2=0.8;                        % Potencia mecanica de la mï¿½quina 2
xg2=0.2;                        % Impedancia del generador 2


%------------------------ Datos de la Carga -------------------------------
H3=15;                         % Constante de inercia del motor de inducción
ka=22.62;%22.62;                     % Constante de amortiguacion por cargas rotantes
kc=1.5;                         % Constante de carga
m=0;                          % Constante de carga
kf=0;%0.1;
Re= 0.00257; 
Xe= 0.01276;
Rm= 50;
Xm= 2.04;
Rr=0.00257;
Xr=0.01276;

Ze= Re+1j*Xe;
Zm= inv(inv(Rm)+inv(1j*Xm));
Zr= Rr+1j*Xr;

%-------------------- Creacion de la Ybus----------------------------------

R=load('Data1.txt');
N=3;                        % Numero de barras   
[Ybus, G, B] = BuiltYbus ( N, R );

Rf=load('Data2.txt');
Nf=4;                        % Numero de barras   
[Ybusf, Gf, Bf] = BuiltYbus ( Nf, Rf );
%------------------------- Ybus Aumentada ---------------------------------
% Ybus aumentanda se le introduce las impedancias de los generadores:

Ybus_aum = Ybus; Ybus_aum(1,1) = inv( 1j * xg1 ) + Ybus(1,1);
Ybus_aum(2,2) = inv( 1j* xg2 ) + Ybus(2,2);
Yb = Ybus_aum;
Ybpf = Yb;



Ybus_aumf = Ybusf; Ybus_aumf(1,1) = inv( 1j * xg1 ) + Ybusf(1,1);
Ybus_aumf(2,2) = inv( 1j* xg2 ) + Ybusf(2,2);
Ybf = Ybus_aumf;


%%%%%%%%%%%%%%%%%%%%%%% FLUJO DE CARGA INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No es Generico para simplificar el Algoritmo, ya que es un caso
% especifico
%-------------- Resolucion del sistema de escuaciones ---------------------
x0=[0 h 0 1];
[x,FVAL,EXITFLAG,OUTPUT,JACOB]=fsolve(@(x)PowerFlow(x,G, B, theta1, V1,...
    V2, Pg2, kc, m, Ze, Zm, Zr,kf),x0,options);
[F,Sd3]=PowerFlow(x,G, B, theta1, V1,...
    V2, Pg2, kc, m, Ze, Zm, Zr,kf);
theta2 = x(1); theta3 = x(3); V3 = x(4); s = x(2); 
[Pg1, Qg1, Qg2, Qg3, DS] = PowerBalance (R, G, B, theta1, theta2, theta3, V1,...
    V2, V3, Pg2, Sd3);
Pref=real(Sd3);

%%%%%%%%%%%%%%%%%%% Calculo de la E' de cada Generador %%%%%%%%%%%%%%%%%%%%
% Calculo de la potencias aparente
S1 = Pg1 + 1j * Qg1;
S2 = Pg2 + 1j * Qg2;

Pm1 = Pg1;
Pm2 = Pg2;
% Calculo de las corrientes
I1=conj(S1/(V1*exp(1j*theta1)));
I2=conj(S2/(V2*exp(1j*theta2)));

% Calculo final:
E_1 = V1 * exp( 1j * theta1 ) + 1j * xg1 * I1;
E_2 = V2 * exp( 1j * theta2 ) + 1j * xg2 * I2;

% Modulos de E'y deltas
ME1 = abs( E_1 );
ME2 = abs( E_2 );

delta1 = angle( E_1 );
w1 = 0;
delta2 = angle( E_2 );
w2 = 0;

%------------------ Corriente Inyectada en la Barra -----------------------
I_iny1 = E_1 / ( 1j * xg1 );
I_iny2 = E_2 / ( 1j * xg2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESTUDIO ESTATICO %%%%%%%%%%%%%%%%%%%%%%%%%%%
s_p=nonzeros(0:h*1e-2:1);
x = [ delta1, delta2, w1, w2, s ];
for k=1:numel(s_p)
    x(end)=s_p(k);
    Iiny = [0;0;0];
    [Pe1, Pe2, Pe3p(k), Qe1, Qe2, V_nodalp(:,k), Q1, Q2, Pejep(k), Pmp(k),...
        Qe3, Pf,I3,Vejep(k) ]= Algoritm (x,ME1, ME2, xg1, xg2, Iiny, Yb, Ze,...
        Zm, Zr, kc, m, kf);    
    Iiny = [0;0;0;0];
    [Pe1, Pe2, Pe3pf(k), Qe1, Qe2, V_nodalpf(:,k), Q1, Q2, Pejepf(k), Pmpf(k),...
        Qe3, Pf,I3,Vejepf(k) ]= Algoritm (x,ME1, ME2, xg1, xg2, Iiny, Ybf, Ze,...
        Zm, Zr, kc, m, kf);
end
figure(1);plot(s_p,Pe3p,s_p,Pmp,s_p,Pejep,'LineWidth',2);grid;title('Curva P3 vs s');xlabel('s');ylabel('P');legend('Pd3','Pm','Peje');
figure(2);plot(Pe3p,abs(V_nodalp(3,:)),Pe3pf,abs(V_nodalpf(3,:)),Pref,abs(V_nodalp(3,:)),'k','LineWidth',2);grid;title('Curva V vs P');xlabel('P');ylabel('V');legend('1','2');
if EXITFLAG~=1
    error  ('Amigo Usuario: Carga mecanica demasiado alta no hay solucion posible');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
% ---------------------- Condiciones inciales ------------------------------
y = [ delta1, delta2, w1, w2, s ];
x = y;
Iiny = [0;0;0];
for k = 1 : numel(t1)

    x0 = y;
    [Pe1o1(k), Pe2o1(k), Pe3o1(k), Qe1o1(k), Qe2o1(k), V_nodal1(:,k), Q1,...
        Q2, Peje_o1(k), Pm_o1(k), Qe3o1(k), Pf_o1(k),I3_o1(k),Vejeo1(k) ]= Algoritm (y, ME1, ME2, xg1, xg2, Iiny, Yb,...
        Ze, Zm, Zr,  kc, m, kf );
       
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( @(x)Ecuacion(x,  y, h, ws, H1, H2,...
            ka, Pm1, Pm2, H3, kc, m, ME1, ME2, xg1, xg2, Iiny, Yb, Pe1o1(k),...
            Pe2o1(k), Peje_o1(k), Pm_o1(k),Ze, Zm, Zr, kf, Pf_o1(k)), x0, options);
    end
    % Solucion
    sol1(:,k) = x';
    y = x;
    [DS] = PowerBalance2 (R, G, B, Pe1o1(k), Pe2o1(k), Pe3o1(k), Q1, Q2,...
        Qe3o1(k), V_nodal1(:,k));
    if DS > 1e-10
        break
    end
end
figure(19)
plot(t(1:numel(Pe1o1)),Pe1o1,t(1:numel(Pe1o1)),Pe2o1,t(1:numel(Pe1o1)),Pe3o1)
grid
legend('1','2','3')
title('P')
figure(20)

plot(t(1:numel(Pe1o1)),Qe1o1,t(1:numel(Pe1o1)),Qe2o1,t(1:numel(Pe1o1)),Qe3o1,t(1:numel(Pe1o1)), Qg3)
legend('1','2','3')
grid
title('Q')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ALGORITMO DE CALCULO DE POTENCIA %%%%%%%%%%%%%%%%%%%%%%


Iiny = [0;0;0;0];

for k = 1 : numel(t2)
    
   x0 = y;
     [Pe1o2(k), Pe2o2(k), Pe3o2(k), Qe1o2(k), Qe2o2(k), V_nodal2(:,k), Q1,...
        Q2, Peje_o2(k), Pm_o2(k),Qe3o2(k), Pf_o2(k),I3_o2(k),Vejeo2(k) ]= Algoritm (y, ME1, ME2,...
        xg1, xg2, Iiny, Ybf, Ze, Zm, Zr,  kc, m, kf );
          
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( @(x)Ecuacion(x,  y, h, ws, H1, H2,...
            ka, Pm1, Pm2, H3, kc, m, ME1, ME2, xg1, xg2, Iiny, Ybf, Pe1o2(k),...
            Pe2o2(k), Peje_o2(k), Pm_o2(k),Ze, Zm, Zr, kf, Pf_o2(k)), x0, options);
    end
    % Solucion
    if x(end)>1
        x(end)=1;
    end
    sol2(:,k) = x';
    y = x;
    [DS] = PowerBalance3 (Rf, Gf, Bf, Pe1o2(k), Pe2o2(k), Pe3o2(k), Q1, Q2,Qe3o2(k), V_nodal2(:,k));
    if sum(abs(DS)) > 1e-10
        break
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
Iiny = [0;0;0];

for k = 1 : numel(t3)
  
   x0 = y;
    [Pe1o3(k), Pe2o3(k), Pe3o3(k), Qe1o3(k), Qe2o3(k), V_nodal3(:,k), Q1,...
        Q2, Peje_o3(k), Pm_o3(k),Qe3o3(k), Pf_o3(k), I3_o3(k),Vejeo3(k) ] = Algoritm (y, ME1, ME2,...
        xg1, xg2, Iiny, Yb, Ze, Zm, Zr,  kc, m, kf  );
       
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( @(x)Ecuacion( x, y, h, ws, H1, H2,...
            ka, Pm1, Pm2, H3, kc, m, ME1, ME2, xg1, xg2, Iiny, Yb, Pe1o3(k),...
            Pe2o3(k), Peje_o3(k), Pm_o3(k),Ze, Zm, Zr, kf, Pf_o3(k)), x0, options);
    end
    % Solucion
    if x(end)>1
        x(end)=1;
    end
   
    sol3(:,k) = x';
    y = x;
    [DS] = PowerBalance2 (R, G, B, Pe1o3(k), Pe2o3(k), Pe3o3(k), Q1, Q2,...
       Qe3o3(k), V_nodal3(:,k));
    if sum(abs(DS)) > 1e-10
        break
    end
    
end



sol = [ sol1, sol2, sol3 ];
delta1 = sol(1,:);
delta2 = sol(2,:);
w1     = sol(3,:);
w2     = sol(4,:);
s    = sol(5,:);

V = [ V_nodal1 V_nodal2(1:3,:) V_nodal3 ];
P1 = [ Pe1o1 Pe1o2  Pe1o3];
P2 = [ Pe2o1 Pe2o2 Pe2o3];
P3 = [ Pe3o1 Pe3o2 Pe3o3];
Q3 = [ Qe3o1 Qe3o2 Qe3o3];
Peje= [Peje_o1 Peje_o2 Peje_o3];
Pm= [Pm_o1 Pm_o2 Pm_o3];
I3= [I3_o1 I3_o2 I3_o3];

Veje= [Vejeo1 Vejeo2 Vejeo3];
    
figure(1000);plot(Q3,abs(V(3,:)),'LineWidth',2);grid;title('Curva Qc vs V3');xlabel('Qc');ylabel('V3');
figure(1);plot(Pe3p,abs(V_nodalp(3,:)),'--',Pe3pf,abs(V_nodalpf(3,:)),'--',Pref*ones(size(Pe3p)),abs(V_nodalp(3,:)),'k',P3,abs(V(3,:)),'LineWidth',2);grid;title('Curva Pc vs V3');xlabel('Pc');ylabel('V3');
figure(111);plot(Peje,abs(Veje),'LineWidth',2);grid;title('Curva Pc vs V3');xlabel('Pc');ylabel('V3');
figure(2);plot(s_p,Pejep,'--',s_p,Pejepf,'--',s,Peje,s_p,Pmp,'k','LineWidth',2); grid;title('Curva Pc vs s');xlabel('Deslizamiento s');ylabel('Pc');legend('Pb','Peje');
figure(3);plot(Pe3p,abs(V_nodalp(3,:)),Pejep,abs(Vejep),'LineWidth',2);grid;title('Curva V vs P');xlabel('P');ylabel('V');
figure(3);plot(s,abs(V(3,:)),s,abs(Veje),'LineWidth',2);grid;title('Curva V3 vs s');xlabel('Deslizamiento s');ylabel('|V3|');
figure(4);plot(t(1:numel(V(3,:))),abs(V(3,:)),'LineWidth',2);grid;title('Curva V3 vs t');xlabel('tiempo [t]');ylabel('V3');
figure(5);plot(t,s,'LineWidth',2);grid;title('Curva s vs t');xlabel('tiempo [t]');ylabel('s');
figure(6);plot(t,P3,'LineWidth',2);grid;title('Curva Pc vs t');xlabel('tiempo [t]');ylabel('Pc');
figure(7);plot(t,delta1-delta2,'LineWidth',2);grid;title('Curva delta1-delta2 vs t');xlabel('tiempo [t]');ylabel('delta');
figure(8);plot(t,w1-w2,'LineWidth',2);grid;title('Curva w1-w2 vs t');xlabel('tiempo [t]');ylabel('w'),
figure(9);plot(t,P1,t,P2,'LineWidth',2);grid;title('Curva Pe1 y Pe2 vs t');xlabel('tiempo [t]');ylabel('Pe'),
figure(10);plot(t,abs(I3),'LineWidth',2);grid;title('Curva I3 vs t');xlabel('tiempo [t]');ylabel('I estator'),
figure(11);plot(t,abs(V(1,:)),t,abs(V(2,:)),'LineWidth',2);grid;title('Curva V1 y V2 vs t');xlabel('tiempo [t]');ylabel('V1');
% figure(12);plot(delta(1:end-1)-angle(V3),Pe_g,'LineWidth',2);grid;title('Curva Pe vs delta-theta3');xlabel('delta-theta3');ylabel('Pe');
% figure(13),plot(t(1:end-1),abs(Ig),'LineWidth',2),grid;title('Curva Ig vs t');xlabel('t');ylabel('Ig')
figure(14);plot(t,delta1,t,delta2,'LineWidth',2);grid;title('Curva delta1 delta2 vs t');xlabel('tiempo [t]');ylabel('delta'),legend('delta1','delta2');
figure(15);plot(t,w1,t,w2,'LineWidth',2);grid;title('Curva w1 w2 vs t');xlabel('tiempo [t]');ylabel('omega'),legend('w1','w2');
figure(16);plot(t,abs(V(1,:)),t,abs(V(2,:)),t,abs(V(3,:)),'LineWidth',2);grid;title('Curva V1, V2 y V3 vs t');xlabel('tiempo [t]');ylabel('V'),legend('V1','V2','V3');
figure(17);plot(t,angle(V(1,:)),t,angle(V(2,:)),t,angle(V(3,:)),'LineWidth',2);grid;title('Curva teta1 , teta2, y teta3  vs t');xlabel('tiempo [t]');ylabel('tetas'),legend('teta1','teta2','teta3');

% figure(6);plot(t,P3,t,Pe_v,'LineWidth',2);grid;title('Curva Pc vs t');xlabel('tiempo [t]');ylabel('Pc');
% Cambiar Pc, por P bornes 
toc