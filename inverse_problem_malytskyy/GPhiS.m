function rez=GPhiS;
global VpKv VsKv k w mue s h KilkistShariv NormalorExp Experement; 
%Дл_ півпростору окремо S хвиля
%Шукаємо Gphi5, Gphi6 : rez(1), rez(2)
%В результаті фомуємо вектор rez із 2-ма компонентами
%Викрористовуэмо Власні значенн_ з  exp(-Bt) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Метод Власних значень та векторів дл_ ШАРУВАТОГО СЕРЕДОВИЩА
% Exp(-Bt(1)-Bt(2)-...-Bt(s))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SumaExponent;
SumaExponent = 1;
for qq=1:s
    Bt_mas(qq) = (( k*k - (w*w)/(VsKv(qq)) )^(1/2))/k;  
    BtR =  h(qq)*k*Bt_mas(qq);
    SumaExponent = SumaExponent * (exp(-BtR));
end;
Dobytok = 1;
for qq=1:s-1
    Dyg(qq) = 1+ ( mue(qq+1)*Bt_mas(qq+1) ) / ( mue(qq)*Bt_mas(qq) );
    Dobytok = Dobytok*Dyg(qq);
end;
ind = s;
Bt_s = (( k*k - (w*w)/(VsKv(ind)) )^(1/2))/k;             %betta(i);
G_s=1+Bt_s*Bt_s;      Mu_s = mue(ind);
DribPiMu = 1/(2*pi*Mu_s);
L1 = (Dobytok*SumaExponent)/(2^s);
Resultat=zeros(2,1);
%%% Znahodimo Gphi1
Resultat(1) = -DribPiMu*L1;
%%% Znahodimo Gphi2
Resultat(2) = -DribPiMu*L1/(2*Bt_s);

rez=Resultat;