%% 
% 
% 
% Using the parameter values form Chol_Input_Tiled, with CL = 0.5

clc
clear
k1 = 3;
k2 = 2;
k3 = 5;
k4 = 4;
k5 = 5;
k6 = 1;
k7 = 4;
k8 = 3; 
k9 = 1;
k10 = 1.2;
k11 = 1;
k12 = 9;
k13 = 1;
k14 = 1;
k15 = 1;
k16 = 1;
p1 = 30;
p2 = 4;
p3 = 1.5;
mu = 25;
eta = 10;
alpha = 10;
theta = 10;
CL = 0.5;

syms Sci C Sr Sh R Cf Cp Ce E HR H Sp P

F=zeros(13,1);
e1=F(1) == (mu - eta*Sci*C);
e2=F(1) == (theta*Ce - eta*Sci*C);
e3=F(3) == (k1*Sci - k2*Sr);
e4=F(4) == (k1*Sci - k10*Sh);
e5=F(5) == (p1*Sr - k11*R - k14*P*R);
e6=F(6) == (k3*CL*R - k4*Cf);
e7=F(7) == (k4*Cf - k5*Cp + k6*Ce);
e8=F(8) == (k5*Cp - k6*Ce + k9*HR*H - k12*Ce - k7*Ce + k8*E);
e9=F(9) == (k7*Ce - k8*E);
e10=F(10) == (p2*Sh - k15*HR*Ce);
e11=F(11) == (alpha - k9*HR*H);
e12=F(12) == (k1*Sci - k13*Sp);
e13=F(13) == (p3*Sp - k16*P);
%% 
% Now solve to find the steady state values in terms of Sci, C, Sr, Sh, R, Cf, 
% Cp, Ce, E, HR, H, Sp, P, Cl and the parameters as listed, then determine the 
% Jacobian J:


[Scieq,Ceq,Sreq,Sheq,Req,Cfeq,Cpeq,Ceeq,Eeq,HReq,Heq,Speq,Peq]=solve([e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13],[Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P]);
f = [rhs(e1),rhs(e2),rhs(e3),rhs(e4),rhs(e5),rhs(e6),rhs(e7),rhs(e8),rhs(e9),rhs(e10),rhs(e11),rhs(e12),rhs(e13)];
J = jacobian(f, [Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P]);
%J = simplify(J)
%% 
% Substitute the steady state values obtained into the Jacobian:


Jeq = subs(J,[Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P],[Scieq,Ceq,Sreq,Sheq,Req,Cfeq,Cpeq,Ceeq,Eeq,HReq,Heq,Speq,Peq])
%% 
% Find the eigenvalue matrix:

double(eig(Jeq))
%% 
% 
% 
% 
% 
%