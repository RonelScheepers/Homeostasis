clear
clc
%Define the symbolic variables and system of ODE for the thirteen species
syms Sci C Sr Sh R Cf Cp Ce E HR H Sp P k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 p1 p2 p3 mu eta alpha theta CL
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

%Now solve to find the steady state values in terms of (Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P,Cl and the parameters as listed), then determine the Jacobian J:
[Scieq,Ceq,Sreq,Sheq,Req,Cfeq,Cpeq,Ceeq,Eeq,HReq,Heq,Speq,Peq]=solve([e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13],[Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P]);
f = [rhs(e1),rhs(e2),rhs(e3),rhs(e4),rhs(e5),rhs(e6),rhs(e7),rhs(e8),rhs(e9),rhs(e10),rhs(e11),rhs(e12),rhs(e13)];
J = jacobian(f, [Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P]);
J = simplify(J) 
%Substitute the steady state values obtained into the Jacobian:
Jeq = subs(J,[Sci,C,Sr,Sh,R,Cf,Cp,Ce,E,HR,H,Sp,P],[Scieq,Ceq,Sreq,Sheq,Req,Cfeq,Cpeq,Ceeq,Eeq,HReq,Heq,Speq,Peq])

% Now run the code for a collection of randomly chosen parameter values (say 10000 or so different parameter sets), and to check how many random
% parameter sets have eigenvalues with a negative real part.

%Number of trials
n = 10000;
found = zeros(n,1);

for i = 1:n
  % Uniformly distributed values for the parameters 
k1num(i) = 1+50*rand();
k2num(i) = 1+50*rand();
k3num(i) = 1+50*rand();
k4num(i) = 1+50*rand();
k5num(i) = 1+50*rand();
k6num(i) = 1+50*rand();
k7num(i) = 1+50*rand();
k8num(i) = 1+50*rand();
k9num(i) = 1+50*rand();
k10num(i) = 1+50*rand();
k11num(i) = 1+50*rand();
k12num(i) = 1+50*rand();
k13num(i) = 1+50*rand();
k14num(i) = 1+50*rand();
k15num(i) = 10+50*rand();
k16num(i) = 1+50*rand();
p1num(i) = 1+50*rand();
p2num(i) = 25+50*rand();
p3num(i) = 1+50*rand();
munum(i) = 1+50*rand();
etanum(i) = 1+50*rand();
alphanum(i) = 1+50*rand();
thetanum(i) = 1+50*rand();
CLnum(i) = 1+50*rand();

 Jeqnum = subs(Jeq,[k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 p1 p2 p3 mu eta alpha theta CL], ...
      [k1num(i) k2num(i) k3num(i) k4num(i) k5num(i) k6num(i) k7num(i) k8num(i) k9num(i) k10num(i) k11num(i) k12num(i) ...
      k13num(i) k14num(i) k15num(i) k16num(i) p1num(i) p2num(i) p3num(i) munum(i) etanum(i) alphanum(i) thetanum(i) CLnum(i)]);

  % Compute eigenvalues
  v{i} = double(eig(Jeqnum));

  % Check if all real parts are < 0
  if any(real(v{i})>=0) 
    found(i) = 1;
  end
end

% Lists the number of trials where eigenvalues with real part >= 0 occurred.
sum(found)
%(found)%The 1 indicates the position of the parameter set (out of 100) that contains one or more eigenvalues >=0.
% The actual eigenvalues can be viewed by opening 'v' in the workspace.