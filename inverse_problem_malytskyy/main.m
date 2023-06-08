function main;
clear all;
% % ������ ��������� �������� )))))))))))))))))))))))))))))))))))))))))))
% ChoisKomponenta - ����� ���_��_ (����������) ������: (P-SV, SH)
%             ChoisKomponenta = 0  - �� ��� ����������: U_z, U_r, U_phi
%             ChoisKomponenta = 1  - ����� ���_��_:     U_z
%             ChoisKomponenta = 2  - ����� ���_��_:     U_r
%             ChoisKomponenta = 3  - ���� ���_��_:     U_phi

% ChoisIntegral -  ������ �������� ������: (P-SV, SH)
%             ChoisIntegral = 0  - ���� ���� ���������: [0, w/Vp]+[w/Vp,w/Vs]
%             ChoisIntegral = 1  - ������ ��������:      [0, w/Vp]
%             ChoisIntegral = 2  - ������ ��������:      [w/Vp, w/Vs]
%             ChoisIntegral = 12 - ��� ��������� ����� �:[0, w/Vs]

% ChiselnijIntegral -  ����� ������������ �� �:
%   !!         ChiselnijIntegral = 1  - ��_��������� (P-SV, SH)
%              ChiselnijIntegral = 2  - ѳ������ (P-SV)
%              ChiselnijIntegral = 3  - ����� � ������ ������� (P-SV)

% KrokIntegr   -  ����� �������_ ��� ����������� �� �: (P-SV, SH)
%             KrokIntegr = 1 - ʲ��ʲ�� ����� �������_ ����� 
%                              ���� �������_ ������� ���-��_ �̲����   
%
%             KrokIntegr = 2 - ʲ��ʲ��� ����� �������_ ������Ӫ���� 
%                              ���� �������_ ������� ���-��_ ������ 

% Metod_Rec_Nad_Pid - ����� ���������_ ������� (P-SV)
%             Metod_Rec_Nad_Pid = 0 - �����������
%             Metod_Rec_Nad_Pid = 1 - ��� ��������
%             Metod_Rec_Nad_Pid = 2 - ϲ� ��������
%   !!!       Metod_Rec_Nad_Pid = 3 - ϲ� ��������, 6�6 � ���������� �� ����������
%             Metod_Rec_Nad_Pid = 4 - ��� ��������, ����� �������_
% MetodPHI_Nad_Pid  - ����� ���������_ ������� ��_ ���������� Uphi (SH)
%                      MetodPHI_Nad_Pid = 1 - ��� �������� !! �� ������ �������
%                      MetodPHI_Nad_Pid = 2 - ϲ� ��������

% IF_4x4_OR_6x6  - ����� ���������_ ����������  (P-SV)  
%             IF_4x4_OR_6x6 = 0   -   4�4 �������
%             IF_4x4_OR_6x6 = 1   -   6�6 �������
% Zona_FIN = =Zona_Far_Intern_Near - � _�� ��� ���������� ����
%             Zona_Far_Near = 0 - �� ����
%             Zona_Far_Near = 1 - �����_   { ���  r}
%             Zona_Far_Near = 2 - ������� { �  1/r} + �����_{�  1/rr}
%POLE      - _�� ���� ������ 
%             POLE=0  - ����� �������� ������� �� ����������� ���������
%             POLE=1  - ������ P �� S  �� ���������� ���������
ChoisKomponenta   = 0;
ChoisIntegral     = 12;
ChiselnijIntegral = 1;
KrokIntegr        = 1;
Metod_Rec_Nad_Pid = 3; %!!��_ ������� �� ���_ ����� ��������������� 6�6
MetodPHI_Nad_Pid  = 1;
IF_4x4_OR_6x6     = 0;
Zona_Far_Near = 1;
POLE = 0
Inverse = 1;
disp(' ---------------- Main Function ------------------ ');
clc;
global  tI Hw UxPfft UxSfft UyPfft UySfft UzPfft UzSfft VpKv VsKv k w mue h s KilkistShariv wI m Vp Vs KilRozbuttjaInt KytPhi r N Fmax t Ht ;
%*********************************************************************************
%***********************  ��������� ����ײ  **************************************
%**********************************************************************************
format long e;
KilkistShariv =3      
s =2 
nFur =11            %�����-���_ ��_ �������_ ������� length(wI) = 2^nFur
KilRozbuttjaInt =2000 %�����-���_ ��_ �����������_ �� �




%*********************************************************************************
%***********************  �������� ������   *************************************
%*********************************************************************************
   
%%_______ ������ ������ ����������  � ������ ������_____________

N=37666

%_________������� ������_ �����������_ ����
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
%_________������ ������_ �� ����


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
% %________________ �� ������� ���� �� ����� � ������ ���______________
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

%����� � ����
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
