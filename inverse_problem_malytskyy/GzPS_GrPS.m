function rez=GzPS_GrPS(ChoisKomponenta,IF_4x4_OR_6x6);
global VpKv VsKv k w mue s h KilkistShariv Experement;
%Дл_ півпростору окремо P та S хвилі
%Шукаємо Gz1P, Gz2P, Gz3P, Gr1P, Gr2P, Gr3P: rez(1)...rez(6)
%        Gz1S, Gz2S, Gz3S, Gr1S, Gr2S, Gr3S: rez(7)...rez(12)
%В результаті фомуємо вектор rez із 12-ма компонентами
%Викрористовуэмо Власні значенн_ з exp(-Al), exp(-Bt) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Метод Власних значень та векторів дл_ ШАРУВАТОГО СЕРЕДОВИЩА\
% Exp(-Al(1)-Al(2)-...-Al(s));  Exp(-Bt(1)-Bt(2)-...-Bt(s))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SumaExponent_Al;  
clear SumaExponent_Bt;
SumaExponent_Al = 1;
SumaExponent_Bt = 1;
for qq=1:s
    Al_mas(qq) = (( k*k - (w*w)/(VpKv(qq)) )^(1/2))/k;  
    AlR =  h(qq)*k*Al_mas(qq);
    SumaExponent_Al = SumaExponent_Al * (exp(-AlR));
    
    Bt_mas(qq) = (( k*k - (w*w)/(VsKv(qq)) )^(1/2))/k;  
    BtR =  h(qq)*k*Bt_mas(qq);
    SumaExponent_Bt = SumaExponent_Bt * (exp(-BtR));
    G_mas(qq) = 1+Bt_mas(qq)*Bt_mas(qq);
end;
Dob_Al = 1;     Dob_Bt = 1;     Drib = 1;
for qq=1:s-1
    kp1 = (mue(qq+1)*G_mas(qq+1)) / (2*mue(qq)) - 1;
    kp2 = (G_mas(qq)) / 2 - (mue(qq+1)) / (mue(qq));
     
    Dyg_Al = kp1 + kp2*( Al_mas(qq+1) )/( Al_mas(qq) ) ;
    Dyg_Bt = kp1 + kp2*( Bt_mas(qq+1) )/( Bt_mas(qq) ) ;
    
    Dob_Al = Dob_Al*Dyg_Al;
    Dob_Bt = Dob_Bt*Dyg_Bt;
    Drib  = Drib/(G_mas(qq)-2);
end;
Drib = Drib/(G_mas(s)-2);
a1 = Drib*Dyg_Al*SumaExponent_Al*Al_mas(1);
a2 = Drib*Dyg_Bt*SumaExponent_Bt*(-1);
b1 = Drib*Dyg_Al*SumaExponent_Al;
b2 = Drib*Dyg_Bt*SumaExponent_Bt*(-Bt_mas(1));


ind = s;
Al_s = (( k*k - (w*w)/(VpKv(ind)) )^(1/2))/k;             %alpha(i);
Bt_s = (( k*k - (w*w)/(VsKv(ind)) )^(1/2))/k;             %betta(i);
G_s=1+Bt_s*Bt_s;      Mu_s = mue(ind);
DribPiMu = 1/(2*pi*Mu_s);

%____________________________________________________________________________
%Шукаємо Gz1P, Gz2P, Gz3P, Gr1P, Gr2P, Gr3P: rez(1)...rez(6)
%%........................  Z  ..........................
%%% Gz1P
    Resultat(1)  = DribPiMu*a1;
%%% Gz2P
    Resultat(2)  =-DribPiMu*a1*Al_s/2;
%%% Gz3P
    Resultat(3)  =DribPiMu*a1/(2*Al_s);
%%........................  R  ..........................
%%% Gr1P
    Resultat(4)  = DribPiMu*b1; 
%%% Gr2P
    Resultat(5) = DribPiMu*b1*Al_s/2;
%%% Gr3P
    Resultat(6)  =- DribPiMu*b1/(2*Al_s);
%____________________________________________________________________________
%Шукаємо Gz1S, Gz2S, Gz3S, Gr1S, Gr2S, Gr3S: rez(7)...rez(12)
%%........................  Z  ..........................
%%% Gz1S
    Resultat(7)  = DribPiMu*a2*G_s/(2*Bt_s);
%%% Gz2S
    Resultat(8)  =-DribPiMu*a2/2; 
%%% Gz3S
    Resultat(9)  = DribPiMu*a2/2;
%%........................  R  ..........................
%%% Gr1S
    Resultat(10)  = DribPiMu*b2*G_s/(2*Bt_s); 
%%% Gr2S
    Resultat(11) = DribPiMu*b2/2; 
%%% Gr3S
    Resultat(12)  =-DribPiMu*b2/2; 

%   КІНЕЦЬ   КІНЕЦЬ    КІНЕЦЬ    КІНЕЦЬ    КІНЕЦЬ    КІНЕЦЬ    КІНЕЦЬ  
rez=Resultat;
