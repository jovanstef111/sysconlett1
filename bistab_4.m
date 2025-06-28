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

G=ss(A,B,C,D);

% Argument transformation of G
bet=-0.05;
B0=B;
A0=(A-bet*eye(n))*(eye(n)-bet*A)^(-1);
C0=(1-bet^2)*C*(eye(n)-bet*A)^(-2);
D0=D+bet*C*(eye(n)-bet*A)^(-1)*B;

% Transformed plant transfer matrix
G0=ss(A0,B0,C0,D0);

 [X,L,F] = care(A0,B0,1*eye(n),1*eye(m)); % needed for coprime fact. of G0
 % Checking
 eig(A0-B0*F) % stable
 % %%%%%%%%%%%%%%%%%%%
 
 % coprime factors
 bQ=ss(A0-B0*F,B0,-F,eye(m));
 bP=ss(A0-B0*F,B0,C0-D0*F,D0);
 
 s=tzero(bP)
 % interpol. points.
 s1=s(1);
 s2=s(9);
 s3=s(10);
 s4=s(7);
 s5=s(8);
 
 ell=5; % number of interpol. points
 
 % interpol. vectors
 [AM,BM,CM,DM]=ssdata(bP);
 x1=null1(DM+CM*(s1*eye(n)-AM)^(-1)*BM);
 x2=null1(DM+CM*(s2*eye(n)-AM)^(-1)*BM);
 x3=null1(DM+CM*(s3*eye(n)-AM)^(-1)*BM);
 x4=null1(DM+CM*(s4*eye(n)-AM)^(-1)*BM);
 x5=null1(DM+CM*(s5*eye(n)-AM)^(-1)*BM);

 [AY,BY,CY,DY]=ssdata(bQ);
 y1=(DY+CY*(s1*eye(n)-AY)^(-1)*BY)*x1;
 y2=(DY+CY*(s2*eye(n)-AY)^(-1)*BY)*x2;
 y3=(DY+CY*(s3*eye(n)-AY)^(-1)*BY)*x3;
 y4=(DY+CY*(s4*eye(n)-AY)^(-1)*BY)*x4; 
 y5=(DY+CY*(s5*eye(n)-AY)^(-1)*BY)*x5; 
 
Api=diag([s1,s2,s3,s4,s5]);
Cmin=[x1,x2,x3,x4,x5];
Cpl=[y1,y2,y3,y4,y5];

% transformation of complex to real interpolat. data
T0=[1 -complex(0,1); 1 complex(0,1)]/sqrt(2);
T=[eye(3) zeros(3,2); 
    zeros(2,3) T0];
hApi=real(T'*Api*T);
hCmin=real(Cmin*T);
hCpl=real(Cpl*T);

% finding the auxiliary Qaux
[Qaux]=nult(hApi,hCmin);

% Iterations to find matrices Q_i
Q(:,:,1)=hCpl;
for i=1:15
    i
[ind,Q(:,:,i+1)]=prv(hApi,Q(:,:,i),hCmin,hCmin); % Qaux or hCmin
    if ind==1
    break
    end
end

% Iterations to find the SPR factors Hm
H=1;
for ii=1:i+1
    Cplii=Q(:,:,ii);
    
    if ii==i+1
    Cminii=hCmin;
    else
    Cminii=Q(:,:,ii+1);
    end
    
XX=lyap(hApi',-(Cminii'*Cplii+Cplii'*Cminii)); % Pick matrix

% checking
e=eig(XX)
% %%%%%%%%

% SPR interpolant
% finding vector xi needed for reduced-order interpolant
      xi1=null((-s1'*eye(ell)+hApi')*XX);
      norm( (Cminii-Cplii)*xi1 ) - norm( (Cminii+Cplii)*xi1 )
      
      xi2=null((-s2'*eye(ell)+hApi')*XX);
      norm( (Cminii-Cplii)*xi2 ) - norm( (Cminii+Cplii)*xi2 )
      
      xi3=null((-s3'*eye(ell)+hApi')*XX);
      norm( (Cminii-Cplii)*xi3 ) - norm( (Cminii+Cplii)*xi3 )
      
      xi4=null((-s4'*eye(ell)+hApi')*XX);
      norm( (Cminii-Cplii)*xi4 ) - norm( (Cminii+Cplii)*xi4 )
      
      xi5=null((-s5'*eye(ell)+hApi')*XX);
      norm( (Cminii-Cplii)*xi5 ) - norm( (Cminii+Cplii)*xi5 )

      if ii~=i+1
      pom1=(Cminii-Cplii)*[xi1,xi3];
      pom2=(Cminii+Cplii)*[xi1,xi3];
      else
      pom1=(Cminii-Cplii)*[xi1,xi3];
      pom2=(Cminii+Cplii)*[xi1,xi3];
      end
      U=real(pom1*pinv(pom2));
      % U=[0.7 0.01; 0.03 0.2]; % checking the formula for interpolant with U
      norm(U)
  
   AH= XX*hApi-((Cplii+Cminii)'+(Cplii-Cminii)'*U)*(eye(m)+U)^(-1)*Cminii;
   BH=((Cplii+Cminii)'+(Cplii-Cminii)'*U)*(eye(m)+U)^(-1);
   CH=Cplii-(eye(m)-U)*(eye(m)+U)^(-1)*Cminii;
   DH=(eye(m)-U)*(eye(m)+U)^(-1);
   EH=XX;
Hm=dss(AH, BH, CH, DH, EH);
Hm=minreal(Hm);

% checking Hm
[AH,BH,CH,DH]=ssdata(Hm);
[eH,eH1]=size(AH);
pole(Hm)
tzero(Hm)
% checking int. condit. on Hm
pom1=Cplii*T';
pom2=Cminii*T';
(DH+CH*(s1*eye(eH)-AH)^(-1)*BH)*pom2(:,1)-pom1(:,1)
(DH+CH*(s2*eye(eH)-AH)^(-1)*BH)*pom2(:,2)-pom1(:,2)
(DH+CH*(s3*eye(eH)-AH)^(-1)*BH)*pom2(:,3)-pom1(:,3)
(DH+CH*(s4*eye(eH)-AH)^(-1)*BH)*pom2(:,4)-pom1(:,4)
(DH+CH*(s5*eye(eH)-AH)^(-1)*BH)*pom2(:,5)-pom1(:,5)

% checking the SPR property of Hm
PP=ispr(AH,BH,CH,DH);
ePP=eig(PP)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=H*Hm; % The final interpolant H is product of i+1 RMs Hm
end

% checking H
[AH,BH,CH,DH,EH]=dssdata(H);
[eH,eH1]=size(AH);
pole(H)
tzero(H)
% checking int. condit on Hm
(DH+CH*(s1*EH-AH)^(-1)*BH)*Cmin(:,1)-Cpl(:,1)
(DH+CH*(s2*EH-AH)^(-1)*BH)*Cmin(:,2)-Cpl(:,2)
(DH+CH*(s3*EH-AH)^(-1)*BH)*Cmin(:,3)-Cpl(:,3)
(DH+CH*(s4*EH-AH)^(-1)*BH)*Cmin(:,4)-Cpl(:,4)
(DH+CH*(s5*EH-AH)^(-1)*BH)*Cmin(:,5)-Cpl(:,5)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cons=(bQ-H)*bP^(-1); % controller (transformed)
Contr=balreal(minreal(minreal(Cons))); % canceling its unstable modes (interp.points)
%Contr=balreal(minreal(minreal(stabsep(Cons)))); %same result with minreal

[AA, BB, CC, DD, EE]=dssdata(Contr);

pp=pole(Contr)  %  Contr. is stable
zz=tzero([ eye(m) G0; Contr eye(m)]) % and CLS is stable

% checking
[pomn,xyz]=size(zz);
pom=zeros(pomn,1);
for iii=1:pomn
pom(iii)=(zz(iii)+bet)/(1+bet*zz(iii));
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inverse argument transformation of the transformed contr.
[AC,BC,CC,DC]=ssdata(Contr);
[nC,nC1]=size(AC);
betm=-bet;
B00=BC;
A00=(AC-betm*eye(nC))*(eye(nC)-betm*AC)^(-1);
C00=(1-betm^2)*CC*(eye(nC)-betm*AC)^(-2);
D00=DC+betm*CC*(eye(nC)-betm*AC)^(-1)*BC;

% Controller
Con=ss(A00,B00,C00,D00); 

% checking the stab. of contr. and of CLS with Con
pp1=pole(Con)
zz1=tzero([ eye(m) G; Con eye(m)])

% reductiom of controller
Conr=reduce(Con,10);

% checking the stab. of contr. and of CLS with Conr
ppr=pole(Conr)  %  kontr. e stabilen
zzr=tzero([ eye(m) G; Conr eye(m)]) % i CLS e stabilen


function [ind,Qopt]=prv(A,Q,Qaux,hCmin)

[m,el]=size(hCmin);

setlmis([])
X=lmivar(1,[el,1]);
X1=lmivar(1,[el,1]);
Y=lmivar(2,[m,el]);
al=lmivar(1,[1,1]);

lmiterm([1 1 1 X],1,A,'s');
lmiterm([1 1 1 Y],Q',-1,'s');

lmiterm([-2 1 1 X],1,1);

lmiterm([3 1 1 X1],1,A,'s');
lmiterm([3 1 1 Y],hCmin',-1,'s');

lmiterm([-4 1 1 X1],1,1);

LMIs=getlmis;

Nn = decnbr(LMIs); 
c = zeros(Nn,1);

for jj=1:Nn, 
	[tauj] = defcx(LMIs,jj,al); 
	c(jj) = tauj; 
end

options=[0,0,0,0,0];
output1 = evalc('[TH,xfeas] = mincx(LMIs,c,options);');
%Xopt=dec2mat(LMIs,xfeas,X); 

% Check feasibility
if isempty(TH)
    ind=0; % disp('No feasible solution found.');
elseif isinf(TH) || isnan(TH)
    ind=0; % disp('Optimization failed or returned an invalid cost.');
else
    ind=1; % disp('Feasible solution found.');
    Qopt=dec2mat(LMIs,xfeas,Y);
    Qopt=Qopt/norm(Qopt);
    return
end


% if the previous code was not succesfull (ind=0):

setlmis([])
X=lmivar(1,[el,1]);
Y=lmivar(2,[m,el]);
al=lmivar(1,[1,1]);

lmiterm([1 1 1 X],1,A,'s');
lmiterm([1 1 1 Y],Q',-1,'s');

lmiterm([-2 1 1 X],1,1);

lmiterm([-3,1,1,al],1,1);
lmiterm([-3,1,2,Y],1,1);
lmiterm([-3,1,2,0],-Qaux);
lmiterm([-3,2,2,al],1,1);

LMIs=getlmis;

Nn = decnbr(LMIs); 
c = zeros(Nn,1);

for jj=1:Nn, 
	[tauj] = defcx(LMIs,jj,al); 
	c(jj) = tauj; 
end

options=[0,0,0,0,0];
output1 = evalc('[TH,xfeas] = mincx(LMIs,c,options);');
Xopt=dec2mat(LMIs,xfeas,X);

Qopt=dec2mat(LMIs,xfeas,Y);

end


function [Qaux]=nult(A,hCmin)

[m,el]=size(hCmin);

setlmis([])
X=lmivar(1,[el,1]);
Y=lmivar(2,[m,el]);
al=lmivar(1,[1,1]);

lmiterm([1 1 1 X],1,A,'s');
lmiterm([1 1 1 Y],hCmin',-1,'s');

lmiterm([-2 1 1 X],1,1);
lmiterm([-2,1,1,al],1,1);

lmiterm([-3,1,1,0],1);
lmiterm([-3,1,2,Y],1,1);
lmiterm([-3,2,2,0],1);

LMIs=getlmis;

Nn = decnbr(LMIs); 
c = zeros(Nn,1);

for jj=1:Nn, 
	[tauj] = defcx(LMIs,jj,al); 
	c(jj) = tauj; 
end

options=[0,0,0,0,0];
output1 = evalc('[TH,xfeas] = mincx(LMIs,c,options);');
Qaux=dec2mat(LMIs,xfeas,Y);

% checking the eig(Xopt)

%Xopt=dec2mat(LMIs,xfeas,X);
%e=eig(Xopt)
%comparison with Qaux=hCmin
%XX=lyap(A',-2*hCmin'*hCmin); % Pick matrix
%eXX=eig(XX)

end

function Popt=ispr(A,B,C,D)

[n,n1]=size(A);

setlmis([])
P=lmivar(1,[n,1]);

lmiterm([1 1 1 P],1,A,'s');
lmiterm([1 1 2 P],1,B);
lmiterm([1 1 2 0],-C');
lmiterm([1 2 2 0],-D-D');
lmiterm([-2 1 1 P],1,1);

LMIs=getlmis;

output1 = evalc('[val,xfeas]=feasp(LMIs);');
Popt=dec2mat(LMIs,xfeas,P);
end