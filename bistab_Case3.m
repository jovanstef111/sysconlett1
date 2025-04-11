% New bistable interpolant and application to strong stabilization
clear

% plant dimensions and matrices
n=10;
m=3;

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
A=A+2*eye(n);
B=-[ -2 -1 0; 1 3 -3; -1 0 2; -2 3 -1; 2 0 1; 
    -1 2 0; -2 -3 -1; 0 1 -2; -1 0 1; -1 3 0];
    
C=-[ -3 0 -1 3 2 1 -2 0 1 -3  ; 
    -2 3 0 -1 -1 -2 -1 -3 1 0 ;
    -2 0 3 1 2 1 -2 -3 0 -1 ];
D=[0 -1 -1; 0 1 2; 0 -3 1];

% plant transfer matrix
G=ss(A,B,C,D);

tic   

% define beta for the transformation
bet=-0.05;
betm=-bet;

% transformation of the plant (A,B,C,D) to (A0,B0,C0,D0) using beta
B0=B;
A0=(A-bet*eye(n))*(eye(n)-bet*A)^(-1);
C0=(1-bet^2)*C*(eye(n)-bet*A)^(-2);
D0=D+bet*C*(eye(n)-bet*A)^(-1)*B;
G0=ss(A0,B0,C0,D0);

% Finding matrices of a coprime factorization of G0 ( =bP*bQ^(-1) )
 [X,L0F,F] = care(A0,B0,1*eye(n),1*eye(m)); 
 bQ=ss(A0-B0*F,B0,-F,eye(m));
 bP=ss(A0-B0*F,B0,C0-D0*F,D0);
 
 ' checking the stability of the coprime factors '
 eig(A0-B0*F)

 ' invariant zeros of the plant '
 s=tzero(bP)

 % Choosing the interpolation points as those s's such that Re(s)\geq 0
 ' interpolation points '
 s1=s(1)
 s2=s(9)
 s3=s(10)
 s4=s(7)
 s5=s(8)
 
 % the number of interpolation conditions/points
 ell=5;
 
 % Creating the interpolation vectors
 [AM,BM,CM,DM]=ssdata(bP);
 x1=null1(DM+CM*(s1*eye(n)-AM)^(-1)*BM);
 x2=null1(DM+CM*(s2*eye(n)-AM)^(-1)*BM);
 x3=null1(DM+CM*(s3*eye(n)-AM)^(-1)*BM);
 x4=null1(DM+CM*(s4*eye(n)-AM)^(-1)*BM);
 x5=conj(x4);

 [AY,BY,CY,DY]=ssdata(bQ);
 y1=(DY+CY*(s1*eye(n)-AY)^(-1)*BY)*x1;
 y2=(DY+CY*(s2*eye(n)-AY)^(-1)*BY)*x2;
 y3=(DY+CY*(s3*eye(n)-AY)^(-1)*BY)*x3;
 y4=(DY+CY*(s4*eye(n)-AY)^(-1)*BY)*x4; 
 y5=(DY+CY*(s5*eye(n)-AY)^(-1)*BY)*x5; 

 % Forming the matrices Api, Cmin and Cpl
Api=diag([s1,s2,s3,s4,s5]);
Cmin=[x1,x2,x3,x4,x5];
Cpl=[y1,y2,y3,y4,y5];

% and transforming them to real ones
T0=[1 -complex(0,1); 1 complex(0,1)]/sqrt(2);
T=[eye(1) zeros(1,2) zeros(1,2); 
    zeros(2,1) eye(2) zeros(2,2);
    zeros(2,1) zeros(2,2) T0];
hApi=real(T'*Api*T);
hCmin=real(Cmin*T);
hCpl=real(Cpl*T);


% Method of Case 3
nH=ell; % order of the interpolant

% Finding matrix Z by minimization of || DH0^(-1)*C10 || 
VN=null(hCpl);
YN=null(VN');
VV=[YN,VN];

pom=hCmin*VV;
pom1=pom(:,1:m);
pom2=pom(:,m+1:ell);
Z=-YN*pom1^(-1);
DHC=[zeros(m,m), pom2]*VV';
pom=[hApi; zeros(nH-ell,ell); -hCpl]*[eye(ell)+Z*hCmin, Z];
pom1=pom(:,1:ell);
pom2=pom(:,ell+1:ell+m);
pomA=[pom1(1:nH,:),zeros(nH,nH-ell)];
pomC=[zeros(nH-ell,ell), eye(nH-ell); hCmin, zeros(m,nH-ell) ];
 
[X,L0,L] = care(pomA',pomC',1*eye(nH),1*eye(m+nH-ell));
L=-L';

 ' checking the stability of H '
 eig(pomA+L*pomC)
 
 V12=L(:,nH-ell+1:nH-ell+m);
 A2=L(:,1:nH-ell);
 
A11=[pom1(1:nH,:)+V12*hCmin, A2];
A21=[pom1(nH+1:nH+m,:),zeros(m,nH-ell)];
A12=pom2(1:nH,:)+V12;
A22=pom2(nH+1:nH+m,:);
B1=zeros(nH,m);
B2=eye(m);
C1=pomC;
C2=[zeros(nH-ell,m); eye(m)];

bA=A11-(A12/A22)*A21;
bB=B1-A12*(A22\B2);
bC=C1-C2*(A22\A21);
bD=-C2*(A22\B2);

% Calling SOF stabilization function
bK=funsof(bA',bC',bB');
bK=bK';

' checking is the SOF stabilization successful '
% i.e. the minimum phase property of H
eig(bA+bB*bK*bC)

K=bK/(eye(m)+bD*bK);

pom=[A11 A12; A21 A22]+[B1; B2]*K*[C1,C2];
AH=pom(1:nH,1:nH);
BH=pom(1:nH,nH+1:nH+m);
CH=pom(nH+1:nH+m,1:nH);
DH=pom(nH+1:nH+m,nH+1:nH+m);

' checking the bistability of H '
eig(AH)
eig(AH-(BH/DH)*CH)

' checking if H satisfies the interpolation conditions '
(DH+CH*(s1*eye(nH)-AH)^(-1)*BH)*x1-y1
(DH+CH*(s2*eye(nH)-AH)^(-1)*BH)*x2-y2
(DH+CH*(s3*eye(nH)-AH)^(-1)*BH)*x3-y3
(DH+CH*(s4*eye(nH)-AH)^(-1)*BH)*x4-y4
(DH+CH*(s5*eye(nH)-AH)^(-1)*BH)*x5-y5

% bistable interpolant
H=ss(AH,BH,CH,DH);

% Finding strong controller for the transformed plant
 Con1=minreal((bQ-H)*bP^(-1));
 
 % Applying on the controller the inverse transformatiion with beta
[AC,BC,CC,DC]=ssdata(Con1);
[nC,nC1]=size(AC);
B00=BC;
A00=(AC-betm*eye(nC))*(eye(nC)-betm*AC)^(-1);
C00=(1-betm^2)*CC*(eye(nC)-betm*AC)^(-2);
D00=DC+betm*CC*(eye(nC)-betm*AC)^(-1)*BC;

% The strong controller for the original plant
Con=ss(A00,B00,C00,D00); 
 
' checking that the controller is stable stabilizing '
pole(Con)
tzero([eye(m) G; Con eye(m)])

% reduction of the controller
Conr=reduce(Con,8);
 
' checking that the reduced controller is stable stabilizing '
pole(Conr)
tzero([eye(m) G; Conr eye(m)])

toc

% function for SOF stabilization by Kimura's method: Find K such that A+B*K*C is stable
function K=funsof(A,B,C)
[n,m]=size(B);
[p,pom]=size(C);
D=null(B');
T=[D,B];
%eig(T);
pomA=T^(-1)*A*T;
pomC=C*T;
A11=pomA(1:n-m,1:n-m);
A12=pomA(1:n-m,n-m+1:n);
A21=pomA(n-m+1:n,1:n-m);
A22=pomA(n-m+1:n,n-m+1:n);
C1=pomC(:,1:n-m);
C2=pomC(:,n-m+1:n);

 [X0,Eig0,L0]=care(A11,A12,1*eye(n-m),1*eye(m));
 L0=-L0;
 %Eig0
 %eig(A11+A12*L0)
 
hA21=-L0*(A11+A12*L0)+A21+A22*L0;
hC1=C1+C2*L0;
K0=-hA21*pinv(hC1);
N=(null(hC1'))';
hC2=N*C2;
hA22=-L0*A12+A22+K0*C2;

[XX,Eig,X]=care(hA22',hC2',0.1*eye(m),5*eye(p+m-n));
X=-X';
%eig(hA22+X*hC2);
 
K=K0+X*N;

end