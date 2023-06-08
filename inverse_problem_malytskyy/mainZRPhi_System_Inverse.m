function rez = mainZRPhi_System_Inverse(ChoisIntegral,ChoisKomponenta,Metod_Rec_Nad_Pid,IF_4x4_OR_6x6, MetodPHI_Nad_Pid, Zona_Far_Near);

global tI Hw UxPfft UxSfft UyPfft UySfft UzPfft UzSfft VpKv VsKv k w mue h s KilkistShariv wI m Vp Vs KilRozbuttjaInt KytPhi r N Fmax t Ht ;
disp(' --------------------- P-SV, SH - випадок ------------------------ ');

    
    Vpj=[3800 4500 4500 6224];
    
    Vsj=[1850 2800 2800 3750];
    Rhoj=[2304 2570 2570 2863];

KytPhij=[5.04];
rj=[2189000];



N1j=[33304];
N2j=[33326];
N3j=[37640];
N4j=[37666];

hj=[10000 12000 2000 24000];  % hs=22 km, s=2, n=3

NN=1
IndVs = 2;%KilkistShariv+1
names1=['Ux_real.dat'];

names2=['Uy_real.dat'];

names3=['Uz_real.dat'];

for jj=1:length(wI)
   w=wI(jj)
for jjj=1:NN %names1(jjj, :)   
UxPP   = load(names1(jjj, :));
UxSS   = load(names1(jjj, :));
UyPP   = load(names2(jjj, :));
UySS   = load(names2(jjj, :));
UzPP   = load(names3(jjj, :));
UzSS   = load(names3(jjj, :));
N1=N1j(jjj, :);
N2=N2j(jjj, :);
N3=N3j(jjj, :);
N4=N4j(jjj, :);

for kj = 1:N1 
    UxPP(kj)= 0;
    UyPP(kj)= 0 ;
    UzPP(kj)= 0;
end
for kj = N1:N2 
    UxPP(kj)= UxPP(kj);
    UyPP(kj)= UyPP(kj);
    UzPP(kj)= UzPP(kj);
end
 for kj = N2:N 
    UxPP(kj)= 0;
    UyPP(kj)= 0 ;
    UzPP(kj)= 0;
 end
 %N=NN
for kjj=1:N
    UxP(kjj)=UxPP(kjj);
    UyP(kjj)=UyPP(kjj);
    UzP(kjj)=UzPP(kjj);
end
for kj = 1:N3 
    UxSS(kj)= 0;
    UySS(kj)= 0 ;
    UzSS(kj)= 0;
end
for kj = N3:N4 
    UxSS(kj)= UxSS(kj);
    UySS(kj)= UySS(kj);
    UzSS(kj)= UzSS(kj);
end
 for kj = N4:N 
    UxSS(kj)= 0;
    UySS(kj)= 0 ;
    UzSS(kj)= 0;
 end
 
for kjj=1:N
    UxS(kjj)=UxSS(kjj);
    UyS(kjj)=UySS(kjj);
    UzS(kjj)=UzSS(kjj);
end

Vp=Vpj(jjj, :);
Vs=Vsj(jjj, :);
Rho=Rhoj(jjj, :);
KytPhi=KytPhij(jjj, :);
r=rj(jjj, :);
h=hj(jjj, :);
VpKv=Vp.^2;
VsKv=Vs.^2;
mue=Rho.*VsKv;

 
    UxPfft1=0;
    UxSfft1=0;
    UyPfft1=0;
    UySfft1=0;
    UzPfft1=0;
    UzSfft1=0;
for qq=1:length(tI)
   t=tI(qq);
   UxPfft1 = UxPfft1+( UxP(qq)*exp(-i*w*t) )*Ht; 
   UxSfft1 = UxSfft1+( UxS(qq)*exp(-i*w*t) )*Ht; 
   UyPfft1 = UyPfft1+( UyP(qq)*exp(-i*w*t) )*Ht; 
   UySfft1 = UySfft1+( UyS(qq)*exp(-i*w*t) )*Ht; 
   UzPfft1 = UzPfft1+( UzP(qq)*exp(-i*w*t) )*Ht; 
   UzSfft1 = UzSfft1+( UzS(qq)*exp(-i*w*t) )*Ht; 
end;%for jj=1:length(wI)   { Кінець Циклу по W }
UxPfft(jj) = UxPfft1;
UxSfft(jj) = UxSfft1;
UyPfft(jj) = UyPfft1;
UySfft(jj) = UySfft1;
UzPfft(jj) = UzPfft1;
UzSfft(jj) = UzSfft1;
%end; %{qq=1:length(tI)}



    U_x_P(jj)  = UxPfft(jj);
    U_y_P(jj)  = UyPfft(jj);
    U_z_P(jj)  = UzPfft(jj);
    
    U_x_S(jj)  = UxSfft(jj);
    U_y_S(jj)  = UySfft(jj);
    U_z_S(jj)  = UzSfft(jj);

%__________________ЦИКЛ ПО W_______________________________________________

   aa = 1/50000000;           cc = w/Vs(IndVs);
   hh = abs(aa-cc)/KilRozbuttjaInt;
%%%________________ІНТЕГРАЛИ по К__________________________________________
   x1=0; x2=0;  x3=0; x4=0;
   x5=0; x6=0;  x7=0; x8=0;
   x9=0; x10=0; x11=0; x12=0;                                                     
%%%________Цикл по К_____________
  for ii=1:KilRozbuttjaInt
        k=aa+(ii-1)*hh;  
        %% --- P - Sv ----  
               GzGr=GzPS_GrPS; %12 штук
               Gz1P = GzGr(1);  Gz2P = GzGr(2);  Gz3P = GzGr(3);
               Gr1P = GzGr(4);  Gr2P = GzGr(5);  Gr3P = GzGr(6);
               Gz1S = GzGr(7);  Gz2S = GzGr(8);  Gz3S = GzGr(9);
               Gr1S = GzGr(10); Gr2S = GzGr(11); Gr3S = GzGr(12);
        %% --- SH     ----
               GPhi=GPhiS;     %2 штуки
               GPhi5 = GPhi(1);  GPhi6 = GPhi(2); 
               z=k*r;
               FunBesselja0=besselj(0,z);
               FunBesselja1=besselj(1,z);
%------------   Інтеграли  дл_ матриці   ----------------------
                x1  = x1+(k*FunBesselja0*Gr1P  * hh);                                               
                x2  = x2+(k*FunBesselja1*Gr2P  * hh);
                x3  = x3+(k*FunBesselja1*Gr3P  * hh);
                x4  = x4+(k*FunBesselja0*Gr1S  * hh);
                x5  = x5+(k*FunBesselja0*GPhi5 * hh);
                x6  = x6+(k*FunBesselja1*Gr2S  * hh);
                x7  = x7+(k*FunBesselja1*GPhi6 * hh);
                x8  = x8+(k*FunBesselja1*Gz1P  * hh);
                x9  = x9+(k*FunBesselja0*Gz2P  * hh);
                x10 = x10+(k*FunBesselja0*Gz3P * hh);
                x11 = x11+(k*FunBesselja1*Gz1S * hh);
                x12 = x12+(k*FunBesselja0*Gz2S * hh);
end;%for ii=1:KilRozbuttjaInt  { Кінець Циклу по К }

   COS1 = cos(KytPhi);
   SIN1 = sin(KytPhi);
   COS2 = cos(2*KytPhi);
   SIN2 = sin(2*KytPhi);
Kmatr(1,1) = x1*COS1*COS1;
Kmatr(1,2) = x1*COS1*SIN1;
Kmatr(1,3) = x2*COS1;
Kmatr(1,4) = x3*COS1*COS1*COS1;
Kmatr(1,5) = x3*COS1*SIN1*SIN1;
Kmatr(1,6) = x3*COS1*SIN2; %!!!!



Kmatr(2,1) = x4*COS1*COS1+x5*SIN1*SIN1;
Kmatr(2,2) =(x4-x5)*COS1*SIN1;
Kmatr(2,3) = x6*COS1;
Kmatr(2,4) =-x6*COS1*COS1*COS1 - x7*SIN1*SIN2;
Kmatr(2,5) =-x6*COS1*SIN1*SIN1 + x7*SIN1*SIN2;
Kmatr(2,6) =-x6*COS1*SIN2 + x7*SIN1*COS2*2;

Kmatr(3,1) = x1*COS1*SIN1;
Kmatr(3,2) = x1*SIN1*SIN1;
Kmatr(3,3) = x2*SIN1;
Kmatr(3,4) = x3*COS1*COS1*SIN1;
Kmatr(3,5) = x3*SIN1*SIN1*SIN1;
Kmatr(3,6) = x3*SIN1*SIN2;


Kmatr(4,1) =(x4-x5)*SIN1*COS1;
Kmatr(4,2) = x4*SIN1*SIN1+x5*COS1*COS1;
Kmatr(4,3) = x6*SIN1;
Kmatr(4,4) =-x6*COS1*COS1*SIN1 + x7*COS1*SIN2;
Kmatr(4,5) =-x6*SIN1*SIN1*SIN1 - x7*COS1*SIN2;
Kmatr(4,6) =-x6*SIN1*SIN2 - x7*COS1*COS2*2;

Kmatr(5,1) = x8*COS1;
Kmatr(5,2) = x8*SIN1;
Kmatr(5,3) = x9;
Kmatr(5,4) = x10*COS1*COS1;
Kmatr(5,5) = x10*SIN1*SIN1;
Kmatr(5,6) = x10*SIN2;


Kmatr(6,1) = x11*COS1;
Kmatr(6,2) = x11*SIN1;
Kmatr(6,3) = x12;
Kmatr(6,4) =-x12*COS1*COS1;
Kmatr(6,5) =-x12*SIN1*SIN1;
Kmatr(6,6) =-x12*SIN2;
%Права частина
Uvect(1) = conj(U_x_P(jj));  
Uvect(2) = conj(U_x_S(jj));   
Uvect(3) = conj(U_y_P(jj));
Uvect(4) = conj(U_y_S(jj));   
Uvect(5) = conj(U_z_P(jj));   
Uvect(6) = conj(U_z_S(jj));



% 

    for l=1:6
    K((jjj-1)*6+l, :)=Kmatr(l, :);
    US((jjj-1)*6+l)=Uvect(l);

end;
end;

KO=K';
U=US';
KK=KO*K;
P=pinv(KK);

MvectREZ=P*KO*U;


%%%_______________Формуємо вектор перетворен Фурє________________________________

    Mxz(jj) = MvectREZ(1);
    Myz(jj) = MvectREZ(2);
    Mzz(jj) = MvectREZ(3);
    Mxx(jj) = MvectREZ(4);
    Myy(jj) = MvectREZ(5);
    Mxy(jj) = MvectREZ(6);
end;%for jj=1:length(wI)   { Кінець Циклу по W }

% %_______________ПЕРЕТВОРЕНН ФУРЄ________________________________
for qq=1:length(tI)
   t = tI(qq)
    furMxz1=0;
    furMyz1=0;
    furMzz1=0;
    furMxx1=0;
    furMyy1=0;
    furMxy1=0;
for jj=1:length(wI)
   w=wI(jj);
furMxz1 = furMxz1+( real(Mxz(jj))*cos(w*t) - imag(Mxz(jj))*sin(w*t) )*Hw;  
furMyz1 = furMyz1+( real(Myz(jj))*cos(w*t) - imag(Myz(jj))*sin(w*t) )*Hw;  
furMzz1 = furMzz1+( real(Mzz(jj))*cos(w*t) - imag(Mzz(jj))*sin(w*t) )*Hw;  
furMxx1 = furMxx1+( real(Mxx(jj))*cos(w*t) - imag(Mxx(jj))*sin(w*t) )*Hw;  
furMyy1 = furMyy1+( real(Myy(jj))*cos(w*t) - imag(Myy(jj))*sin(w*t) )*Hw;  
furMxy1 = furMxy1+( real(Mxy(jj))*cos(w*t) - imag(Mxy(jj))*sin(w*t) )*Hw;  
end;%for jj=1:length(wI)   { Кінець Циклу по W }
    furMxz(qq)=furMxz1/pi;
    furMyz(qq)=furMyz1/pi;
    furMzz(qq)=furMzz1/pi;
    furMxx(qq)=furMxx1/pi;
    furMyy(qq)=furMyy1/pi;
    furMxy(qq)=furMxy1/pi;
end; %{qq=1:length(tI)}
%_____Формуємо результат щоб перекидати в головну програму_______
MM = length(tI)
for jj=1:MM
    U(jj)     =furMxz(jj);
    U(MM+jj)  =furMyz(jj);
    U(2*MM+jj)=furMzz(jj);
    U(3*MM+jj)=furMxx(jj);
    U(4*MM+jj)=furMyy(jj);
    U(5*MM+jj)=furMxy(jj);
end;
rez=U;