function [Peje, Pm, Pf, Sd3,Pd3,Qd3,Id3,Zeq,Veje]=MI(Zr, Zm, Ze,V3, theta3, s,kc,m,kf)

% Motor de Inducción:
Rr= real(Zr);
Reje= Rr*(1-s)/s;
Zrotor= Zr + Reje;
Zeq2= inv(inv(Zm)+inv(Zrotor));
Zeq= Ze + Zeq2;
Ve=V3*exp(1j*theta3);

Id3=Ve/Zeq;
Sd3= abs(Ve)^2/conj(Zeq);
Pd3= real(Sd3);
Qd3= imag(Sd3);


In= Ve/Ze;

Zeq3= inv(inv(Ze)+inv(Zm));

Ieje= In*Zeq3/(Zeq3+Zrotor);

Peje=abs(Ieje)^2*Reje;
Veje= Ieje*Reje;

wr=(1-s);
Pf=kf*wr;
Pm=kc*wr^m;
