% Using hifoo to find a strong controller for the plant of bistab.m
clear

% necessary files for hihoo
addpath('C:\Users\jovans\Desktop\hanso2_0');
addpath('C:\Users\jovans\Desktop\HIFOO3.501');

A=-[2 -3 -2 1 -1 0 1 2 -1 3 ; 
   0 0 -2 1 -1 -2 3 1 1 0  ;
   -2 1 3 2 -3 0 -1 1 3 3  ;
   -3 1 0 -2 0 0 -1 2 3 1  ;
    1 0 -1 0 1 -3 -1 0 0 -1  ;
    -2 3 2 0 -1 0 1 2 -1 -3  ; 
    2 -1 -3 -2 -3 1 0 -1 3 -3  ;
    3 -1 0 -2 0 1 -1 2 3 1  ; 
    2 0 1 3 -1 3 -1 1 0 -1 ;
    -1 0 1 0 -1 3 -1 1 0 -1 ] ;
A=A+2*eye(10);
B=-[ -2 -1 0; 1 3 -3; -1 0 2; -2 3 -1; 2 0 1; 
    -1 2 0; -2 -3 -1; 0 1 -2; -1 0 1; -1 3 0];
C=-[ -3 0 -1 3 2 1 -2 0 1 -3  ; 
    -2 3 0 -1 -1 -2 -1 -3 1 0 ;
    -2 0 3 1 2 1 -2 -3 0 -1 ];
D=[0 -1 -1; 0 1 2; 0 -3 1];

G=ss(A,B,C,D);

% preparation for hifoo
PP.A=A;
PP.B=B;
PP.C=C;
PP.D=D;

Psys={PP,'K'};

ORDER=10; 
INIT=[];
FUN='ss';
UPPERBND= [0,0];
OPTIONS.nrand=100;

% application of the hifoo
[K,val,viol]=hifoo(Psys,ORDER,INIT,FUN,UPPERBND,OPTIONS);

AK=K.a;
BK=K.b;
CK=K.c;
DK=K.d;

KHIFOO=ss(AK,BK,CK,DK);

% checking the poles of the controller and of the closed loo[p system
pole(KHIFOO)
tzero( [eye(3) KHIFOO; G eye(3)])

