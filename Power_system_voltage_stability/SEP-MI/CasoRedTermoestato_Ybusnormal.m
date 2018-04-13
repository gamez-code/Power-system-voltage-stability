%%%%%%%%%%%%%% ALGORITMO DINAMICO DE ESTABILIDAD DE VOLTAJE %%%%%%%%%%%%%%%
%Prof. Juan Bermudez
%Br. Manuel Jimenez
clc
clear all
close all
tic

options=optimset('Display','off');
% Topologia del sistema JB pag 151 con una modificacion, la carga
% no es una potencia constante sino un reostato.

%-------------------- Creacion de la Ybus----------------------------------

%Data:  i  f  R     x    Bshunt
    R=[
        1  2  0.01  0.2  0.1;
        1  3  0.02  0.15 0.2;
        2  3  0.03  0.1  0.3;
        3  3  0     -1.5 0];

% La impedancia de falla, la red equivalente sec - y sec 0 es igual a:    
Zf=inv(0.4555-3.4239*1i);

%Ybus
N=3;                        % Numero de barras
Ybus=zeros(N);              %Se crea una matriz de ceros
nexos=size(R,1);            %el numeros de nexos

for i=1:nexos
   m=R(i,1);%punto de incio
   n=R(i,2);%punto final
       if m-n==0
           %Se introduce la compensacion
           Ybus(m,m)=Ybus(m,m)+inv(R(i,3)+1j*R(i,4))+1j*R(i,5)/2;
       else
           Ybus(m,m)=Ybus(m,m)+inv(R(i,3)+1j*R(i,4))+1j*R(i,5)/2;
           Ybus(n,n)=Ybus(n,n)+inv(R(i,3)+1j*R(i,4))+1j*R(i,5)/2;
           Ybus(m,n)=-inv(R(i,3)+1j*R(i,4));
           Ybus(n,m)=Ybus(m,n);
       end
end           
              % G y B para el flujo de carga
                G=real(Ybus); B=imag(Ybus);


%%%%%%%%%%%%%%%%%%%%%%%% DATOS DEL PROBLEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------Datos del tiempos:-----------------------------------
f = 60;                           % frecuencia
ws = 2 * pi * f;                      % frecuencia angular
T = 1 / f;                          % Periodo
t_i = 100 * T;                      % tiempo de inicio de la falla
t_d = 40 * T;                     % tiempo de despeje
t_f = 500 * T;                     % tiempo final de interacion

h = T * 1e-0;                   % Delta t, pasos de integracion
t1 =   0 : h : t_i;             % tiempo pre-falla
t2 = t_i : h : t_i + t_d;       % tiempo de falla
t3 = t_i + t_d : h : t_i + t_d + t_f;% tiempo post-falla
t = [ t1 t2 t3 ];               % Vector de tiempo


%--------------------Datos de los generadores------------------------------
%Generador 1:
% Tipo Slack:
V1 = 1;                           % Voltaje en bornes 
theta1 = 0;                       % Angulo en bornes
H1 = 100;                         % Constante de inercia de la máquina 1
xg1 = 0.1;                      % Impedancia del generador 1

%Generador 2:
V2 = 1.01;                        % Voltaje en bornes
H2 = 100;                         % Constate de inercia se la máquina 2
Pg2 = 0.8;                        % Potencia mecanica de la máquina 2
xg2 = 0.2;                      % Impedancia del generador 2


%------------------------ Datos de la Carga -------------------------------
tao = 0.10;                      % Constate de tiempo del Reostato
ka = 0.06;
Pref = 3.6;

%%%%%%%%%%%%%%%%%%%%%%% FLUJO DE CARGA INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------- Incognitas ----------------------------------

theta2 = @(x) x(1);
theta3 = @(x) x(2);
V3 = @(x) x(3);
Gc = @(x) x(4);

%---------------------------- Ecuaciones ----------------------------------
%                     ------  Ecuacion 1: ------
%cosenos y senos:
c23 = @(x) cos( theta2(x) - theta3(x) );
s23 = @(x) sin( theta2(x) - theta3(x) );
c21 = @(x) cos( theta2(x) - theta1 );
s21 =@(x) sin( theta2(x) - theta1 );

G22V2 = G(2,2) * V2^2;

P23 = @(x) V2 * V3(x) * ( G(2,3) * c23(x) + B(2,3) * s23(x) );
P21 = @(x)  V2 * V1 * ( G(2,1) * c21(x) + B(2,1) * s21(x) );
SP2k = @(x) P21(x) + P23(x);

Pn2 = @(x) G22V2 + SP2k(x);
DP2 =@(x) Pg2 - Pn2(x);

f1 = @(x) DP2(x);
%                   ------ Ecuacion 2: ------
%cosenos y senos:
c32=@(x) cos(theta3(x)-theta2(x));
s32=@(x) sin(theta3(x)-theta2(x));
c31=@(x) cos(theta3(x)-theta1);
s31=@(x) sin(theta3(x)-theta1);

G33V3=@(x)G(3,3)*V3(x)^2;

P32= @(x)V3(x)*V2*(G(3,2)*c32(x)+B(3,2)*s32(x));
P31= @(x)V3(x)*V1*(G(3,1)*c31(x)+B(3,1)*s31(x));
SP3k=@(x) P31(x)+P32(x);

Pn3=@(x) G33V3(x)+SP3k(x);
DP3=@(x) -Pref-Pn3(x);

f2=@(x) DP3(x);

%                   ------ Ecuacion 3: ------
%cosenos y senos:

B33V3=@(x) B(3,3)*V3(x)^2;

Q32=@(x) V3(x)*V2*(G(3,2)*s32(x)-B(3,2)*c32(x));
Q31=@(x) V3(x)*V1*(G(3,1)*s31(x)-B(3,1)*c31(x));
SQ3k=@(x) Q31(x)+Q32(x);

Qn3=@(x)-B33V3(x)+SQ3k(x);
DQ3=@(x)-Qn3(x);

f3=@(x) DQ3(x);

%                   ------ Ecuacion 4: ------
% Conductancia:
f4=@(x) -Pref+V3(x)^2*Gc(x);

%----------------- Sistemas de ecuaciones a resolver ----------------------
F=@(x)[f1(x) f2(x) f3(x) f4(x)];

%-------------- Resolucion del sistema de escuaciones ---------------------
x0=[0 0 1 0];
[x,FVAL,EXITFLAG,OUTPUT,JACOB]=fsolve(F,x0,options);
% Con esta nueva modificacion usando las funciones handles se optimizo el
% codigo, haciendo que se ejecute en 0.7 segundo, cuando antes lo hacia en
% 3.5 seg
theta2 = x(1); theta3 = x(2); V3 = x(3); Gca = x(4); 

%%%%%%%%%%%%%%%%%%% Calculo de la E' de cada Generador %%%%%%%%%%%%%%%%%%%%
% Calculo de la Potencia activa
c23= cos(theta2-theta3);
s23= sin(theta2-theta3);

c13= cos(theta1-theta3);
s13= sin(theta1-theta3);

c21= cos(theta2-theta1);
s21= sin(theta2-theta1);

c12= cos(theta1-theta2);
s12= sin(theta1-theta2);

G11V1 = G(1,1) * V1^2;
G22V2 = G(2,2) * V2^2;

P13 = V1 * V3 * ( G(1,3) * c13 + B(1,3) * s13 );
P12 = V1 * V2 * ( G(1,2) * c12 + B(1,2) * s12 );

P23 = V2 * V3 * ( G(2,3) * c23 + B(2,3) * s23 );
P21 = V2 * V1 * ( G(2,1) * c21 + B(2,1) * s21 );

SP1k = P12 + P13;
SP2k = P21 + P23;

P1 = G11V1 + SP1k;
P2 = G22V2 + SP2k;


% Calculo de la Potencia reactiva
B11V1 = B(1,1) * V1^2;
B22V2 = B(2,2) * V2^2;

Q23 = V2 * V3 * ( G(2,3) * s23 - B(2,3) * c23 );
Q13 = V1 * V3 * ( G(1,3) * s13 - B(1,3) * c13 );

Q21 = V2 * V1 * ( G(2,1) * s21 - B(2,1) * c21 );
Q12 = V1 * V2 * ( G(1,2) * s12 - B(1,2) * c12 );

SQ1k= Q13 + Q12;
SQ2k= Q23 + Q21;

Q1= -B11V1 + SQ1k;
Q2= -B22V2 + SQ2k;

% Calculo de la potencias aparente
S1 = P1 + 1j * Q1;
S2 = P2 + 1j * Q2;

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


% Para la Ybus de falla le agregamos la impedancia de falla
Ybusf = Ybus_aum; Ybusf(1,1) = Ybusf(1,1) + inv( Zf );
Yb = Ybus_aum;


%%%%%%%%%%%%% PROCESO ITERATIVO PARA CALCULO DE ESTABILIDAD %%%%%%%%%%%%%%%
%---------------------- Condiciones inciales ------------------------------
y = [ delta1, delta2, w1, w2, Gca ];
x = y;
%--------------------------- Incognitas -----------------------------------
delt1 = @(x) x(1);
delt2 = @(x) x(2);
omeg1 = @(x) x(3);
omeg2 = @(x) x(4);
Gc    = @(x) x(5);

%%%%%%%%%%%%%%%%%%% ALGORITMO DE CALCULO DE POTENCIA %%%%%%%%%%%%%%%%%%%%%%
E1 = @(y) ME1 * exp ( 1j * delt1(y) );
E2 = @(y) ME2 * exp ( 1j * delt2(y) );

I_iny1 = @(y) E1(y) / ( 1j * xg1 );
I_iny2 = @(y) E2(y) / ( 1j * xg2 );

Iiny = @(y) [ I_iny1(y); I_iny2(y); 0 ];

Y = @(y,Yb) [Yb(1,1)  Yb(1,2) Yb(1,3); 
             Yb(2,1)  Yb(2,2) Yb(2,3);
             Yb(3,1)  Yb(3,2) Yb(3,3)+Gc(y)];
         
V_nodal = @(y,Yb) inv( Y(y,Yb) ) * Iiny(y);
V1 = @(y,Yb) nonzeros( V_nodal(y,Yb) .* [ 1; 0; 0] );
V2 = @(y,Yb) nonzeros( V_nodal(y,Yb) .* [ 0; 1; 0] );
V3 = @(y,Yb) nonzeros( V_nodal(y,Yb) .* [ 0; 0; 1] );

Vxg1 = @(y,Yb) E1(y) - V1(y,Yb);
Vxg2 = @(y,Yb) E2(y) - V2(y,Yb);

Ig1 = @(y,Yb) Vxg1(y,Yb) / ( 1j * xg1 );
Ig2 = @(y,Yb) Vxg2(y,Yb) / ( 1j * xg2 );

Pe1 = @(y,Yb) real( E1(y) * conj( Ig1(y,Yb) ) );
Pe2 = @(y,Yb) real( E2(y) * conj( Ig2(y,Yb) ) );
Qe1 = @(y,Yb) imag( E1(y) * conj( Ig1(y,Yb) ) );
Qe2 = @(y,Yb) imag( E2(y) * conj( Ig2(y,Yb) ) );

Pe3 = @(y,Yb) abs( V3(y,Yb) )^2 * Gc(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1 : numel(t1)

    F1 = @(x) - omeg1(x) + omeg1(y) + h * ws / ( 4 * H1 ) * ...
        ( P1 - Pe1(x,Yb) - ka * omeg1(x) + P1 - Pe1(y,Yb) - ka * omeg1(y) );
    F2 = @(x) - omeg2(x) + omeg2(y) + h * ws / ( 4 * H2 ) * ...
        ( P2 - Pe2(x,Yb) - ka * omeg2(x) + P2 - Pe2(y,Yb) - ka * omeg2(y) );
    
    F3 = @(x) - delt1(x) + delt1(y) + h / 2 *( omeg1(x) + omeg1(y) );
    F4 = @(x) - delt2(x) + delt2(y) + h / 2 *( omeg2(x) + omeg2(y) );
    
    F5 = @(x) - Gc(x) + Gc(y) + h / ( 2 * tao ) * ( Pref - Pe3(x,Yb) + Pref - Pe3(y,Yb) );
    
    % Sistemas de ecuaciones
    SE = @(x) [ F1(x), F2(x), F3(x), F4(x), F5(x)];
    
    % Condiciones iniciales
    x0 = y;
    % Proceso de resolucion
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( SE, x0, options);
    end
    % Solucion
    sol1(:,k) = x';
    y = x;
    V_1(:,k) = V_nodal(x,Yb);
    P_11(k) = Pe1(x,Yb);
    P_21(k) = Pe2(x,Yb);
    Q_11(k) = Qe1(x,Yb);
    Q_21(k) = Qe2(x,Yb);
    P_31(k) = Pe3(x,Yb);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Yb = Ybusf;

for k = 1 : numel(t2)
    
    F1 = @(x) - omeg1(x) + omeg1(y) + h * ws / ( 4 * H1 ) * ...
        ( P1 - Pe1(x,Yb) - ka * omeg1(x) + P1 - Pe1(y,Yb) - ka * omeg1(y) );
    F2 = @(x) - omeg2(x) + omeg2(y) + h * ws / ( 4 * H2 ) * ...
        ( P2 - Pe2(x,Yb) - ka * omeg2(x) + P2 - Pe2(y,Yb) - ka * omeg2(y) );
    
    F3 = @(x) - delt1(x) + delt1(y) + h / 2 *( omeg1(x) + omeg1(y) );
    F4 = @(x) - delt2(x) + delt2(y) + h / 2 *( omeg2(x) + omeg2(y) );
    
    F5 = @(x) - Gc(x) + Gc(y) + h / ( 2 * tao ) * ( Pref - Pe3(x,Yb) + Pref - Pe3(y,Yb) );
    
    % Sistemas de ecuaciones
    SE = @(x) [ F1(x), F2(x), F3(x), F4(x), F5(x)];


    
    % Condiciones iniciales
    x0 = y;
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( SE, x0, options);
    end
    % Solucion
    sol2(:,k) = x';
    y = x;
    V_2(:,k) = V_nodal(x,Yb);
    
    P_12(k) = Pe1(x,Yb);
    P_22(k) = Pe2(x,Yb);
    Q_12(k) = Qe1(x,Yb);
    Q_22(k) = Qe2(x,Yb);
    P_32(k) = Pe3(x,Yb);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POST-FALLA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yb = Ybus_aum;


for k = 1 : numel(t3)
  
  
    F1 = @(x) - omeg1(x) + omeg1(y) + h * ws / ( 4 * H1 ) * ...
        ( P1 - Pe1(x,Yb) - ka * omeg1(x) + P1 - Pe1(y,Yb) - ka * omeg1(y) );
    F2 = @(x) - omeg2(x) + omeg2(y) + h * ws / ( 4 * H2 ) * ...
        ( P2 - Pe2(x,Yb) - ka * omeg2(x) + P2 - Pe2(y,Yb) - ka * omeg2(y) );
    
    F3 = @(x) - delt1(x) + delt1(y) + h / 2 *( omeg1(x) + omeg1(y) );
    F4 = @(x) - delt2(x) + delt2(y) + h / 2 *( omeg2(x) + omeg2(y) );
    
    F5 = @(x) - Gc(x) + Gc(y) + h / ( 2 * tao ) * ( Pref - Pe3(x,Yb) + Pref - Pe3(y,Yb) );
    
    % Sistemas de ecuaciones
    SE = @(x) [ F1(x), F2(x), F3(x), F4(x), F5(x)];


    
    % Condiciones iniciales
    x0 = y;
    if k == 1
        x = x;
    else
        [ x, FVAL, EXITFLAG ] = fsolve( SE, x0, options);
    end
    % Solucion
    sol3(:,k) = x';
    y = x;
    V_3(:,k) = V_nodal(x,Yb);
    
    P_13(k) = Pe1(x,Yb);
    P_23(k) = Pe2(x,Yb);
    Q_13(k) = Qe1(x,Yb);
    Q_23(k) = Qe2(x,Yb);
    P_33(k) = Pe3(x,Yb);
    
end



sol = [ sol1, sol2, sol3 ];
delta1 = sol(1,:);
delta2 = sol(2,:);
w1     = sol(3,:);
w2     = sol(4,:);
Gca    = sol(5,:);

V = [ V_1 V_2 V_3 ];
P1 = [ P_11 P_12 P_13];
P2 = [ P_21 P_22 P_23];
P3 = [ P_31 P_32 P_33];

figure(1);plot(P3,abs(V(3,:)),'LineWidth',2);grid;title('Curva Pc vs V3');xlabel('Pc');ylabel('V3');
figure(2);plot(Gca,P3,'LineWidth',2); grid;title('Curva Pc vs G');xlabel('Conductancia G');ylabel('Pc');
figure(3);plot(Gca,abs(V(3,:)),'LineWidth',2);grid;title('Curva V3 vs G');xlabel('Conductancia G');ylabel('|V3|');
figure(4);plot(t,abs(V(3,:)),'LineWidth',2);grid;title('Curva V3 vs t');xlabel('tiempo [t]');ylabel('V3');
figure(5);plot(t,Gca,'LineWidth',2);grid;title('Curva G vs t');xlabel('tiempo [t]');ylabel('G');
figure(6);plot(t,P3,'LineWidth',2);grid;title('Curva Pc vs t');xlabel('tiempo [t]');ylabel('Pc');
figure(7);plot(t,delta1-delta2,'LineWidth',2);grid;title('Curva delta1-delta2 vs t');xlabel('tiempo [t]');ylabel('delta');
figure(8);plot(t,w1-w2,'LineWidth',2);grid;title('Curva w1-w2 vs t');xlabel('tiempo [t]');ylabel('w'),
figure(9);plot(t,P1,t,P2,'LineWidth',2);grid;title('Curva Pe1 y Pe2 vs t');xlabel('tiempo [t]');ylabel('Pe'),
% figure(10);plot(t,Q,'LineWidth',2);grid;title('Curva Qe vs t');xlabel('tiempo [t]');ylabel('Qe'),
figure(11);plot(t,abs(V(1,:)),t,abs(V(2,:)),'LineWidth',2);grid;title('Curva V1 y V2 vs t');xlabel('tiempo [t]');ylabel('V1');
% figure(12);plot(delta(1:end-1)-angle(V3),Pe_g,'LineWidth',2);grid;title('Curva Pe vs delta-theta3');xlabel('delta-theta3');ylabel('Pe');
% figure(13),plot(t(1:end-1),abs(Ig),'LineWidth',2),grid;title('Curva Ig vs t');xlabel('t');ylabel('Ig')

toc