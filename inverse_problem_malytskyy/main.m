function main;
clear all;
% % Задаємо ПАРАМЕТРИ ПРОГРАМИ )))))))))))))))))))))))))))))))))))))))))))
% ChoisKomponenta - котре рівн_нн_ (КОМПОНЕНТУ) рахуємо: (P-SV, SH)
%             ChoisKomponenta = 0  - всі три компоненти: U_z, U_r, U_phi
%             ChoisKomponenta = 1  - перше рівн_нн_:     U_z
%             ChoisKomponenta = 2  - друге рівн_нн_:     U_r
%             ChoisKomponenta = 3  - третє рівн_нн_:     U_phi

% ChoisIntegral -  котрий ІНТЕГРАЛ рахуємо: (P-SV, SH)
%             ChoisIntegral = 0  - сума двох інтегралів: [0, w/Vp]+[w/Vp,w/Vs]
%             ChoisIntegral = 1  - перший інтеграл:      [0, w/Vp]
%             ChoisIntegral = 2  - другий інтеграл:      [w/Vp, w/Vs]
%             ChoisIntegral = 12 - два інтеграла разом у:[0, w/Vs]

% ChiselnijIntegral -  метод ІНТЕГРУВАННЯ по К:
%   !!         ChiselnijIntegral = 1  - Пр_мокутників (P-SV, SH)
%              ChiselnijIntegral = 2  - Сімпсона (P-SV)
%              ChiselnijIntegral = 3  - Гауса з трьома точками (P-SV)

% KrokIntegr   -  спосіб РОЗБИТТ_ при інтегруванні по К: (P-SV, SH)
%             KrokIntegr = 1 - КІЛЬКІСТ точок розбитт_ СТАЛА 
%                              КРОК розбитт_ проміжку інт-нн_ ЗМІННИЙ   
%
%             KrokIntegr = 2 - КІЛЬКІСТЬ точок розбитт_ ЗБІЛЬШУЄТЬСЯ 
%                              КРОК розбитт_ проміжку інт-нн_ СТАЛИЙ 

% Metod_Rec_Nad_Pid - метод обчисленн_ МАТРИЦЬ (P-SV)
%             Metod_Rec_Nad_Pid = 0 - Рекурентний
%             Metod_Rec_Nad_Pid = 1 - НАД джерелом
%             Metod_Rec_Nad_Pid = 2 - ПІД джерелом
%   !!!       Metod_Rec_Nad_Pid = 3 - ПІД джерелом, 6х6 в чисельнику та знаменнику
%             Metod_Rec_Nad_Pid = 4 - НАД джерелом, власні значенн_
% MetodPHI_Nad_Pid  - метод обчисленн_ МАТРИЦЬ дл_ компоненти Uphi (SH)
%                      MetodPHI_Nad_Pid = 1 - НАД джерелом !! Це кращий випадок
%                      MetodPHI_Nad_Pid = 2 - ПІД джерелом

% IF_4x4_OR_6x6  - метод обчисленн_ ЗНАМЕННИКА  (P-SV)  
%             IF_4x4_OR_6x6 = 0   -   4х4 матриці
%             IF_4x4_OR_6x6 = 1   -   6х6 матриці
% Zona_FIN = =Zona_Far_Intern_Near - в _кій зоні обчислюємо поле
%             Zona_Far_Near = 0 - всі зони
%             Zona_Far_Near = 1 - дальн_   { БЕЗ  r}
%             Zona_Far_Near = 2 - проміжна { З  1/r} + ближн_{З  1/rr}
%POLE      - _ке поле рахуємо 
%             POLE=0  - повна хвильова картина за нормальними формулами
%             POLE=1  - шматки P та S  за виведеними формулами
ChoisKomponenta   = 0;
ChoisIntegral     = 12;
ChiselnijIntegral = 1;
KrokIntegr        = 1;
Metod_Rec_Nad_Pid = 3; %!!Дл_ повного хв пол_ треба використовувати 6х6
MetodPHI_Nad_Pid  = 1;
IF_4x4_OR_6x6     = 0;
Zona_Far_Near = 1;
POLE = 0
Inverse = 1;
disp(' ---------------- Main Function ------------------ ');
clc;
global  tI Hw UxPfft UxSfft UyPfft UySfft UzPfft UzSfft VpKv VsKv k w mue h s KilkistShariv wI m Vp Vs KilRozbuttjaInt KytPhi r N Fmax t Ht ;
%*********************************************************************************
%***********************  ПАРАМЕТРИ ЗАДАЧІ  **************************************
%**********************************************************************************
format long e;
KilkistShariv =3      
s =2 
nFur =11            %Викор-тьс_ дл_ розбитт_ частоти length(wI) = 2^nFur
KilRozbuttjaInt =2000 %Викор-тьс_ дл_ інтегруванн_ по К




%*********************************************************************************
%***********************  ОБЕРНЕНА ЗАДАЧА   *************************************
%*********************************************************************************
   
%%_______ БЕРЕМО ШМАТКИ СЕЙСМОГРАМ  в часовій області_____________

N=37666

%_________Задаэмо розбит_ перетворенн_ Фурэ
Ht=0.05;
 Fmax1=0.8; %0.3; %2.4;%10.0;%1.;%10.;%1.0;%2.0;%10.0;%10.0;
Fmin1=0.2;
Tmax=N*Ht
Fmax=1/(2*Ht)
 Hf=1/(2*Tmax)
% Hf=1/Tmax


Hw=2*pi*Hf;
Wmax=2*pi*Fmax;
Wmax1=2*pi*Fmax1
Wmin1=2*pi*Fmin1;

 wI=[Wmin1:Hw:Wmax1];
%_________Задаємо розбит_ по часу


tI=[Ht:Ht:Tmax];


Tmax
%_______________________-----------------------------
TenM= mainZRPhi_System_Inverse(ChoisIntegral,ChoisKomponenta,Metod_Rec_Nad_Pid,IF_4x4_OR_6x6, MetodPHI_Nad_Pid, Zona_Far_Near);
    for jj=1:N
Mxz(jj)=TenM(jj);
Myz(jj)=TenM(N+jj) ;
Mzz(jj)=TenM(2*N+jj); 
Mxx(jj)=TenM(3*N+jj);
Myy(jj)=TenM(4*N+jj);
Mxy(jj)=TenM(5*N+jj); 
    end;
% %________________ Всі малюнки один під іншим в одному вікні______________
Mal_Pidrjad=figure;
figure(Mal_Pidrjad);
%%----------
subplot(6,1,1);
hold;   grid;   title( 'M_xz ');
plot (tI,Mxz,'r');
%%----------
subplot(6,1,2);
hold;   grid;   title( 'M_yz ');
plot (tI,Myz,'g');
%%----------
subplot(6,1,3);
hold;   grid;   title( 'M_zz ');
plot (tI,Mzz,'b');
%%----------
subplot(6,1,4);
hold;   grid;   title( 'M_xx ');
plot (tI,Mxx,'y');
%%----------
subplot(6,1,5);
hold;   grid;   title( 'M_yy ');
plot (tI,Myy,'k');
%%----------
subplot(6,1,6);
hold;   grid;   title( 'M_xy ');
plot (tI,Mxy,'m');

%Запис у файл
 r2 = 't, s \n';
 r3 = 'M(t)*, N-m \n';
 r4 = '04 05 2022     23 00 00.0\n';
 r5 = '0 0 \n';
 %r6 = '11167 5.0000000000000E-0001 6 \n';
 r6 = '37666  5.000000000000000000E-0002 6 \n';
% r2 = '';
% r3 = '';
% r4 = '';
% r5 = '';
% r6 = '';

[FMxz,mes]   = fopen('Mxz.arr', 'w');
    r1 = 'Eigenvalues, ~~~ Mxz ~~~ component\n';
    fprintf(FMxz,r1);fprintf(FMxz,r2);fprintf(FMxz,r3);
    fprintf(FMxz,r4);fprintf(FMxz,r5);fprintf(FMxz,r6);
[FMyz,mes]   = fopen('Myz.arr', 'w');
    r1 = 'Eigenvalues, ~~~ Myz ~~~ component \n';
    fprintf(FMyz,r1);fprintf(FMyz,r2);fprintf(FMyz,r3);
    fprintf(FMyz,r4);fprintf(FMyz,r5);fprintf(FMyz,r6);
[FMzz,mes]   = fopen('Mzz.arr', 'w');
    r1 = 'Eigenvalues, ~~~ Mzz ~~~ component \n';
    fprintf(FMzz,r1);fprintf(FMzz,r2);fprintf(FMzz,r3);
    fprintf(FMzz,r4);fprintf(FMzz,r5);fprintf(FMzz,r6);
[FMxx,mes]   = fopen('Mxx.arr', 'w');
    r1 = 'Eigenvalues, ~~~ Mxx ~~~ component \n';
    fprintf(FMxx,r1);fprintf(FMxx,r2);fprintf(FMxx,r3);
    fprintf(FMxx,r4);fprintf(FMxx,r5);fprintf(FMxx,r6);
[FMyy,mes]   = fopen('Myy.arr', 'w');
    r1 = 'Eigenvalues, ~~~ Myy ~~~ component\n';
    fprintf(FMyy,r1);fprintf(FMyy,r2);fprintf(FMyy,r3);
    fprintf(FMyy,r4);fprintf(FMyy,r5);fprintf(FMyy,r6);
[FMxy,mes]   = fopen('Mxy.arr', 'w');
    r1 = 'Eigenvalues, ~~~ Mxy ~~~ component \n';
    %r1 = '';
    fprintf(FMxy,r1);fprintf(FMxy,r2);fprintf(FMxy,r3);
    fprintf(FMxy,r4);fprintf(FMxy,r5);fprintf(FMxy,r6);
for jj=1:N
    %Mxz(200)
    fprintf(FMxz,'%16.16e', Mxz(jj));
    fprintf(FMxz,'\n');
    fprintf(FMyz,'%16.16e', Myz(jj));
    fprintf(FMyz,'\n');
    fprintf(FMzz,'%16.16e', Mzz(jj));
    fprintf(FMzz,'\n');
    fprintf(FMxx,'%16.16e', Mxx(jj));
    fprintf(FMxx,'\n');
    fprintf(FMyy,'%16.16e', Myy(jj));
    fprintf(FMyy,'\n');
    fprintf(FMxy,'%16.16e', Mxy(jj));
    fprintf(FMxy,'\n');
end;

fclose(FMxz);
fclose(FMyz);
fclose(FMzz);
fclose(FMxx);
fclose(FMyy);
fclose(FMxy);
%end;
% 
% 
% 
