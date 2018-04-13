%%%%%%%%%%%%%% ALGORITMO DINAMICO DE ESTABILIDAD DE VOLTAJE %%%%%%%%%%%%%%%
%Prof. Juan Bermudez
%Br. Manuel Jimenez
clc
clear all
close all
tic

options=optimset('Display','off');

%Este es un programa donde se hace un estudio de un GS-Linea-Termostato


%%%%%%%%%%%%%%%%%%%%%%%% DATOS DEL PROBLEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Datos del tiempos:-----------------------------------
f=60;                           % frecuencia
ws=2*pi*f;                      % frecuencia angular
T=1/f;                          % Periodo
t_i=120*T;                       % tiempo de inicio de la falla
t_d=91*T;   %tiempo critico 50 ciclos estabilidad de voltaje y angulo                    % tiempo de despeje
%Nota tiempo critico de despeje con solo estabilidad de angulo es de 38
%ciclos, insensible a variaciones del reostato. Constate de tiempo del
%reostato igual a infinito.
t_f=1000*T;                      % tiempo final de interacion

h=T;                            % Delta t, pasos de integracion
t1=0:h:t_i;                     % tiempo pre-falla
t2=t_i:h:t_i+t_d;               % tiempo de falla
t3=t_i+t_d:h:t_i+t_d+t_f;       % tiempo post-falla
t=[t1 t2 t3];                   % Vector de tiempo


%--------------------Datos de los generadores------------------------------
%Generador 1:
% Tipo Slack:
V1=1;                           % Voltaje en bornes 
theta1=0;                       % Angulo en bornes
H1=10;                          % Constante de inercia de la m�quina 1
xg1=0.1;                        % Impedancia del generador 1
xg1n=0.1;
xg1o=0.05;

%------------------------ Datos de la Carga -------------------------------
tau=0.2;                           % Constante tiempo
Pref=1.5;                          % Potencia de referencia.
ka=0.06*377;

%-------------------- Creacion de la Ybus----------------------------------

R=load('Data1.txt');
N=2;                        % Numero de barras   
[Ybus,G,B]=BuiltYbus(N,R);

% Ybus de falla
Rf=load('Data2.txt');
N=3;                        % Numero de barras   
[Ybusf,Gf,Bf]=BuiltYbus(N,Rf);

% Ybus de Sec. Neg.
Rn=load('Data3.txt');
N=3;                        % Numero de barras   
[Ybusn,Gn,Bn]=BuiltYbus(N,Rn);

% Ybus de Sec. 0
Ro=load('Data4.txt');
N=3;                        % Numero de barras   
[Ybuso,Go,Bo]=BuiltYbus(N,Ro);


%%%%%%%%%%%%%%%%%%%%%%% FLUJO DE CARGA INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No es Generico para simplificar el Algoritmo, ya que es un caso
% especifico
%-------------- Resolucion del sistema de escuaciones ---------------------
x0=[0 1 0];

[x,FVAL,EXITFLAG,OUTPUT,JACOB]=fsolve(@(x)PowerFlow(x, G, B, theta1, V1,...
    Pref),x0,options);
V2=x(2);
theta2=x(1);
Gt=x(3);
[Pg1, Qg1, DS] = PowerBalance (R, G, B, theta1, theta2, V1, V2,Pref);


%%%%%%%%%%%%%%%%%%% Calculo de la E' de cada Generador %%%%%%%%%%%%%%%%%%%%
% Calculo de la potencias aparente
S1 = Pg1 + 1j * Qg1;

Pm1 = Pg1;

% Calculo de las corrientes
I1=conj(S1/(V1*exp(1j*theta1)));

% Calculo final:
E_1 = V1 * exp( 1j * theta1 ) + 1j * xg1 * I1;


% Modulos de E'y deltas
ME1 = abs( E_1 );


delta1 = angle( E_1 );
w1 = 0;


%------------------ Corriente Inyectada en la Barra -----------------------
I_iny1 = E_1 / ( 1j * xg1 );


%------------------------- Ybus Aumentada ---------------------------------
% Ybus aumentanda se le introduce las impedancias de los generadores:

Ybus_aum = Ybus; Ybus_aum(1,1) = inv( 1j * xg1 ) + Ybus(1,1);


% perturbacion
Ybus_aumf = Ybusf; Ybus_aumf(1,1) = inv( 1j * xg1 ) + Ybusf(1,1);
Ybus_aumn = Ybusn; Ybus_aumn(1,1) = inv( 1j * xg1n ) + Ybusn(1,1);
Ybus_aumo = Ybuso; Ybus_aumo(1,1) = inv( 1j * xg1o ) + Ybuso(1,1);


Yb = Ybus_aum;
Ybf= Ybus_aumf;
Ybn= Ybus_aumn;
Ybo= Ybus_aumo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESTUDIO ESTATICO %%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1;
for k=0:0.1:20
    Gtp(m)=k;
    x(3)=Gtp(m);
    Iiny = [0;0];
    [Ybb]=BuiltYbf(Yb,0,0,Gtp(m));
    [Pe1,Ptp(m), Qe1, V1, V2p(m)]= Algoritm (x,ME1, xg1, Iiny, Ybb);
    Iiny = [0;0;0];
    [Ybff]=BuiltYbf(Ybf,Ybn,Ybo,Gtp(m));
    [Pe1,Ptfp(m), Qe1, V1, V2fp(m)]= Algoritm (x,ME1, xg1, Iiny, Ybff);
     m=m+1;
end

Pref_v=Pref*ones(size(V2p));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
% ---------------------- Condiciones inciales ------------------------------
y = [ delta1, w1, Gt ];
x = y;
Iiny = [0;0];
[sol1,Pe1o1,Pto1, V1_1, V2_1,Ig1_1,E1_1]=DinamicSolver(t1,y,x, h, ws, H1,ka,...
            Pm1, ME1, xg1, Iiny, Yb, Pref,tau,options,0,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ALGORITMO DE CALCULO DE POTENCIA %%%%%%%%%%%%%%%%%%%%%%

% Falla monofasica en la mitad de la linea 1-2
%-------------------- Creacion de la Ybus de Falla ------------------------
y = sol1(:,end);
x = y;
Iiny = [0;0;0];
[sol2,Pe1o2,Pto2, V1_2, V2_2,Ig1_2,E1_2,V]=DinamicSolver(t2,y,x, h, ws, H1,ka,...
            Pm1, ME1, xg1, Iiny, Ybf, Pref,tau,options,Ybn,Ybo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
y = sol2(:,end);
x = y;
Iiny = [0;0];
[sol3,Pe1o3,Pto3, V1_3, V2_3,Ig1_3,E1_3]=DinamicSolver(t3,y,x, h, ws, H1,ka,...
            Pm1, ME1, xg1, Iiny, Yb, Pref,tau,options,0,0);



sol = [ sol1, sol2 ,sol3];
delta1 = sol(1,:);
w1     = sol(2,:);
Gt    = sol(3,:);

V2 = [ V2_1 V2_2 V2_3 ];
V1 = [ V1_1 V1_2 V1_3 ];
P1 = [ Pe1o1 Pe1o2 Pe1o3];
Pt = [ Pto1 Pto2 Pto3];
Ig1= [Ig1_1 Ig1_2 Ig1_3];
E1=[E1_1 E1_2 E1_3];

grid minor
figure(1);plot(Ptp,abs(V2p),'--',Ptfp,abs(V2fp),'--',Pref_v,abs(V2p),'k',Pt,abs(V2),'LineWidth',2);title('Curva P vs V2');xlabel('P');ylabel('V2');legend('1','2','3','4')
grid on
grid minor
figure(2);plot(Gtp,Ptp,'--',Gtp,Ptfp,'--',Gtp,Pref_v,'k',Gt,Pt,'LineWidth',2);title('Curva P-G');xlabel('G');ylabel('P');legend('1','2','3','4')
grid on
grid minor
figure(3);plot(t(1:numel(V2)),abs(V2),'LineWidth',2);title('Curva V2 vs t');xlabel('tiempo [t]');ylabel('V2');
grid on
grid minor
figure(4);plot(t,Gt,'LineWidth',2);title('Curva G vs t');xlabel('tiempo [t]');ylabel('G');
grid on
grid minor
figure(5);plot(t,Pt,'LineWidth',2);title('Curva P2 vs t');xlabel('tiempo [t]');ylabel('P2');
grid on
grid minor
figure(6);plot(t,delta1,'LineWidth',2);title('Curva delta1 vs t');xlabel('tiempo [t]');ylabel('delta');
grid on
grid minor
figure(7);plot(t,w1,'LineWidth',2);title('Curva w1 vs t');xlabel('tiempo [t]');ylabel('w'),
figure(8);plot(t,P1,t,Pt,'LineWidth',2);title('Curva P1 vs t');xlabel('tiempo [t]');ylabel('P1');
legend('1','2')
grid on
grid minor
figure(9);plot(t,abs(Ig1),'LineWidth',2);title('Curva Ig1 vs t');xlabel('tiempo [t]');ylabel('Ig1');
grid on
grid minor
figure(10);plot(t(1:numel(V1)),abs(V1),'LineWidth',2);title('Curva V1 vs t');xlabel('tiempo [t]');ylabel('V1');
grid on
grid minor

toc