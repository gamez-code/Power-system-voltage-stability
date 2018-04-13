%%%%%%%%%%%%%% ALGORITMO DINAMICO DE ESTABILIDAD DE VOLTAJE %%%%%%%%%%%%%%%
%Prof. Juan Bermudez
%Br. Manuel Jimenez
clc
clear all
close all
tic

options=optimset('Display','off');

% Use la topologia del sistema de los apuntes JB de
% la Pag 151, con una modificaci�n, la carga no es 
% una potencia constante sino un reostato.
% 
% Corr� el flujo de carga, Hice una verificaci�n 
% haciendo una balance de potencia.


%-------------------- Creacion de la Ybus----------------------------------

%Data:  i  f  R     x    Bshunt
R=load('Data1.txt');
Rf=load('Data2.txt');
N=3;                        % Numero de barras   
Nf=4;
[Ybus, G, B] = BuiltYbus ( N, R );
[Ybusf,Gf,Bf]=BuiltYbus(Nf,Rf);
%%%%%%%%%%%%%%%%%%%%%%%% DATOS DEL PROBLEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Datos del tiempos:-----------------------------------
f=60;                           % frecuencia
ws=2*pi*f;                      % frecuencia angular
T=1/f;                          % Periodo
t_i=10*T;                       % tiempo de inicio de la falla
t_d=49*T;                       % tiempo de despeje
%Nota tiempo critico de despeje con solo estabilidad de angulo es de 38
%ciclos, insensible a variaciones del reostato. Constate de tiempo del
%reostato igual a infinito.
t_f=2880*T;                      % tiempo final de interacion

h=T;                       % Delta t, pasos de integracion
t1=0:h:t_i;                     % tiempo pre-falla
t2=t_i:h:t_i+t_d;               % tiempo de falla
t3=t_i+t_d:h:t_i+t_d+t_f;       % tiempo post-falla
t=[ t1 t2 t3 ];                 % Vector de tiempo


%--------------------Datos de los generadores------------------------------
%Generador 1:
% Tipo Slack:
V1=1;                           % Voltaje en bornes 
theta1=0;                       % Angulo en bornes
H1=10;                          % Constante de inercia de la m�quina 1
xg1=0.1;                        % Impedancia del generador 1

%Generador 2:
% Tipo P-V
V2=1.01;                        % Voltaje en bornes
H2=5;                           % Constate de inercia se la m�quina 2
Pg2=0.8;                        % Potencia mecanica de la m�quina 2
xg2=0.2;                        % Impedancia del generador 2


%------------------------ Datos de la Carga -------------------------------
tau=.2;%0.1;                     % Constate de tiempo del Reostato
% Buscar constantes de tiempo de termoestato
ka=22.62;                      % Constante de amortiguacion por cargas rotantes
Pref=3;                     % Potencia de referencia del Reostato

%%%%%%%%%%%%%%%%%%%%%%% FLUJO DE CARGA INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No es Generico para simplificar el Algoritmo, ya que es un caso
% especifico
%-------------- Resolucion del sistema de escuaciones ---------------------
x0=[0 0 1 0];
[x,FVAL,EXITFLAG,OUTPUT,JACOB]=fsolve(@(x)PowerFlow(x,G, B, theta1, V1, V2, Pg2, Pref),x0,options);
theta2 = x(1); theta3 = x(2); V3 = x(3); Gca = x(4); 
[Pg1, Qg1, Qg2, DS] = PowerBalance (R, G, B, theta1, theta2, theta3, V1, V2, V3, Pg2, Pref);

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

%------------------------- Ybus Aumentada ---------------------------------
% Ybus aumentanda se le introduce las impedancias de los generadores:

Ybus_aum = Ybus; Ybus_aum(1,1) = inv( 1j * xg1 ) + Ybus(1,1);
Ybus_aum(2,2) = inv( 1j* xg2 ) + Ybus(2,2);
Yb = Ybus_aum;


Ybus_aumf = Ybusf; Ybus_aumf(1,1) = inv( 1j * xg1 ) + Ybusf(1,1);
Ybus_aumf(2,2) = inv( 1j* xg2 ) + Ybusf(2,2);
Ybf = Ybus_aumf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESTUDIO ESTATICO %%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1;
x = [ delta1, delta2, w1, w2, Gca ];
for k=0:0.1:250
    Gtp(m)=k;
    x(5)=Gtp(m);
    Iiny = [0;0;0];
     [Pe1, Pe2, Pe3p(m), Qe1, Qe2, V_nodalp(:,m), ...
         I3 ] = Algoritm (x, ME1, ME2, xg1, xg2, Iiny, Yb );
    Iiny = [0;0;0;0];
     [Pe1, Pe2, Pe3pf(m), Qe1, Qe2, V_nodalpf(:,m), ...
         I3 ] = Algoritm (x, ME1, ME2, xg1, xg2, Iiny, Ybf );
    m=m+1;
end

Pref_v=Pref*ones(size(V_nodalp(3,:)));
% plot(Pe3p,abs(V_nodalp(3,:)),Pe3pf,abs(V_nodalpf(3,:)),Pref_v,abs(V_nodalp(3,:)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
% ---------------------- Condiciones inciales ------------------------------
y = [ delta1, delta2, w1, w2, Gca ];
x = y;
Iiny = [0;0;0];
[sol1,Pe1o1, Pe2o1, Pe3o1, Qe1o1, Qe2o1, V_nodal1, ...
         I3o1,I2o1,I1o1]=DinamicSolver(t1,y,x, h, ws, H1,H2,ka,...
            Pm1,Pm2, ME1, ME2, xg1, xg2, Iiny, Yb, Pref,tau,options,R,G,B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ALGORITMO DE CALCULO DE POTENCIA %%%%%%%%%%%%%%%%%%%%%%

y=sol1(:,end);
x=y;
Iiny = [0;0;0;0];
[sol2,Pe1o2, Pe2o2, Pe3o2, Qe1o2, Qe2o2, V_nodal2, ...
         I3o2,I2o2,I1o2]=DinamicSolver(t2,y,x, h, ws, H1,H2,ka,...
            Pm1,Pm2, ME1, ME2, xg1, xg2, Iiny, Ybf, Pref,tau,options,Rf,Gf,Bf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%

Iiny = [0;0;0];
y=sol2(:,end);
x=y;
[sol3,Pe1o3, Pe2o3, Pe3o3, Qe1o3, Qe2o3, V_nodal3, ...
         I3o3,I2o3,I1o3]=DinamicSolver(t3,y,x, h, ws, H1,H2,ka,...
            Pm1,Pm2, ME1, ME2, xg1, xg2, Iiny, Yb, Pref,tau,options,R,G,B);


sol = [ sol1, sol2, sol3 ];
delta1 = sol(1,:);
delta2 = sol(2,:);
w1     = sol(3,:);
w2     = sol(4,:);
Gca    = sol(5,:);


V = [ V_nodal1 V_nodal2(1:3,:) V_nodal3 ];
P1 = [ Pe1o1 Pe1o2 Pe1o3];
P2 = [ Pe2o1 Pe2o2 Pe2o3];
P3 = [ Pe3o1 Pe3o2 Pe3o3];
I3 = [ I3o1 I3o2 I3o3];
I2 = [ I2o1 I2o2 I2o3];
I1 = [ I1o1 I1o2 I1o3];

ngtp=find(Gtp==50);
figure(1);plot(Pe3p,abs(V_nodalp(3,:)),'--',Pe3pf,abs(V_nodalpf(3,:)),'--',Pref_v,abs(V_nodalp(3,:)),'k',P3,abs(V(3,:)),'LineWidth',2);grid;title('Curva P3 vs V3');xlabel('P3');ylabel('V3');
figure(2);plot(Gtp(1:ngtp),Pe3p(1:ngtp),'--',Gtp(1:ngtp),Pe3pf(1:ngtp),'--',Gtp(1:ngtp),Pref_v(1:ngtp),'k',Gca,P3,'LineWidth',2); grid;title('Curva P3 vs G');xlabel('Conductancia G');ylabel('P3');
% figure(3);plot(Gca,abs(V(3,:)),'LineWidth',2);grid;title('Curva V3 vs G');xlabel('Conductancia G');ylabel('|V3|');
figure(4);plot(t,abs(V(3,:)),'LineWidth',2);grid;title('Curva V3 vs t');xlabel('tiempo [t]');ylabel('V3');
% figure(5);plot(t,Gca,'LineWidth',2);grid;title('Curva G vs t');xlabel('tiempo [t]');ylabel('G');
figure(6);plot(t,P3,'LineWidth',2);grid;title('Curva Pc vs t');xlabel('tiempo [t]');ylabel('Pc');
figure(7);plot(t,delta1-delta2,'LineWidth',2);grid;title('Curva delta1-delta2 vs t');xlabel('tiempo [t]');ylabel('delta');
figure(8);plot(t,w1-w2,'LineWidth',2);grid;title('Curva w1-w2 vs t');xlabel('tiempo [t]');ylabel('w'),
figure(9);plot(t,P3,t,P1,t,P2,'LineWidth',2);grid;title('Curva Pe1 y Pe2 vs t');xlabel('tiempo [t]');ylabel('Pe'),
legend('P3','P1','P2')
% figure(10);plot(t,I3,'LineWidth',2);grid;title('Curva I3 vs t');xlabel('tiempo [t]');ylabel('I3'),
% figure(11);plot(t,abs(V(1,:)),t,abs(V(2,:)),'LineWidth',2);grid;title('Curva V1 y V2 vs t');xlabel('tiempo [t]');ylabel('V1');
% figure(12);plot(delta(1:end-1)-angle(V3),Pe_g,'LineWidth',2);grid;title('Curva Pe vs delta-theta3');xlabel('delta-theta3');ylabel('Pe');
% figure(13),plot(t(1:end-1),abs(Ig),'LineWidth',2),grid;title('Curva Ig vs t');xlabel('t');ylabel('Ig')
figure(14);plot(t,delta1,t,delta2,'LineWidth',2);grid;title('Curva delta1 delta2 vs t');xlabel('tiempo [t]');ylabel('delta'),legend('delta1','delta2');
figure(15);plot(t,w1,t,w2,'LineWidth',2);grid;title('Curva w1 w2 vs t');xlabel('tiempo [t]');ylabel('omega'),legend('w1','w2');
figure(99);plot(t,abs(I3),t,abs(I2),t,abs(I1),'LineWidth',2);grid;title('Curva Ig vs t');xlabel('tiempo [t]');ylabel('I');legend('3','2','1');
% figure(16);plot(t,abs(V(1,:)),t,abs(V(2,:)),t,abs(V(3,:)),'LineWidth',2);grid;title('Curva V1, V2 y V3 vs t');xlabel('tiempo [t]');ylabel('V'),legend('V1','V2','V3');
% figure(17);plot(t,angle(V(1,:)),t,angle(V(2,:)),t,angle(V(3,:)),'LineWidth',2);grid;title('Curva teta1 , teta2, y teta3  vs t');xlabel('tiempo [t]');ylabel('tetas'),legend('teta1','teta2','teta3');

% % figure(6);plot(t,P3,t,Pe_v,'LineWidth',2);grid;title('Curva Pc vs t');xlabel('tiempo [t]');ylabel('Pc');
% 
toc