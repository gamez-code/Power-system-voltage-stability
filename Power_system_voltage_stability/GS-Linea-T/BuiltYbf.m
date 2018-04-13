function [Ybf]=BuiltYbf(Ybf,Ybn,Ybo,Gt)
if abs(det(Ybn))==0 && abs(det(Ybo))==0
     Ybf(2,2)=Ybf(2,2)+Gt;
else
    Ybn(2,2)=Ybn(2,2)+Gt;
    Ybo(2,2)=Ybo(2,2)+Gt;
    Ybf(2,2)=Ybf(2,2)+Gt;
    
    Zbo=inv(Ybo);
    Zbn=inv(Ybn);
    
    Ybf(3,3)=Ybf(3,3)+inv(Zbn(3,3)+Zbo(3,3));
end
end
