function [Ybus, G, B] = BuiltYbus (NodesNumber,Branch )
N = NodesNumber;
R = Branch;
%Ybus

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
