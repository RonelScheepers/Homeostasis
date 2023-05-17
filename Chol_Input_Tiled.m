
clear
clc

% Set the initial values
u(1) = 5;
u(2) = 1;
u(3) = 1;
u(4) = 1;
u(5) = 1;
u(6) = 50;
u(7) = 1;
u(8) = 100;
u(9) = 5;
u(10) = 1;
u(11) = 1;
u(12) = 1;
u(13) = 1;

%% Define anonymous function input1(t)
%  Time: 0<=t<50 50<=t<100  100<=t<150  150<=t<200  200<=t  
%  Input:  0.5      1          1.5          1         0.5
input1 = @(t) 0.5 + 0.5*(50<=t & t<100) + 1*(100<=t & t<150) + 0.5*(150<=t & t<200);

% Set the model parameters
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
k16 = 30;
k17 = 4;
k18 = 1.5;
k19 = 25;
k20 = 10;
k21 = 10;
k22 = 10;
k23 = 1;

%% Perform the numerical integration
k = [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23];
init = [u(1) u(2) u(3) u(4) u(5) u(6) u(7) u(8) u(9) u(10) u(11) u(12) u(13)];
tspan = [0 250];
[t,u] = ode23s(@(t,u) gene(t,u,k,input1), tspan, init);

%% Plot results
t1=0:.5:250;
tiledlayout(4,1)

nexttile 
    hold on
plot(t1,input1(t1),'-k', 'LineWidth',2.0)
title(''); legend('CL(Input)')
xlabel('Time t'); ylabel('Concentration')
hold off
ylim([0 2])



nexttile 
    hold on
plot(t,u(:,6),'-',t,u(:,7),'-',t,u(:,8),'-',t,u(:,9),'-', 'LineWidth',2.0)
title(''); legend('Cf','Cp','Ce','E')
xlabel('Time t'); ylabel('Concentration')
ylim([0 6])


nexttile 
    hold on
plot(t,u(:,2),'-',t,u(:,5),'-',t,u(:,11),'-', 'LineWidth',2.0)
title(''); legend('C','R','H')
xlabel('Time t'); ylabel('Concentration')
ylim([0 60])



nexttile 
    hold on
plot(t,u(:,1),'-',t,u(:,3),'-',t,u(:,10),'-',t,u(:,12),'-',t,u(:,13),'-',t,u(:,4),'-', 'LineWidth',2.0)
title('','FontSize',16); legend('Sci','Sr','HR','Sp','P','Sh')
xlabel('Time t'); ylabel('Concentration')
ylim([0 1.5])



%% ODE System
function eqns = gene(t,u,k,in1)
% Using u = [Sci C Sr Sh R Cf Cp Ce  E  HR  H  Sp  P] 
% Equation    1  2  3  4 5  6  7  8  9  10  11 12  13

%Using k17=p1, k18=p2,k19=p3, k20=mu, k21=eta, k22=alpha, k23=theta

eqns = zeros(13,1); % To start with we have 13 empty equations
eqns(1) = k(19) - k(20)*u(1)*u(2);
eqns(2) = k(22)*u(8) - k(20)*u(1)*u(2);
eqns(3) = k(1)*u(1) - k(2)*u(3);
eqns(4) = k(1)*u(1) - k(10)*u(4);
eqns(5) = k(16)*u(3) - k(11)*u(5) - k(14)*u(13)*u(5);
eqns(6) = k(3)*u(5)*in1(t) - k(4)*u(6);
eqns(7) = k(4)*u(6) - k(5)*u(7) + k(6)*u(8);
eqns(8) = k(5)*u(7) - k(6)*u(8) + k(9)*u(10)*u(11) - k(12)*u(8) - k(7)*u(8) + k(8)*u(9);
eqns(9) = k(7)*u(8) - k(8)*u(9);
eqns(10) = k(17)*u(4) - k(15)*u(10)*u(8);
eqns(11) = k(21) - k(9)*u(10)*u(11);
eqns(12) = k(1)*u(1) - k(13)*u(12);
eqns(13) = k(18)*u(12) - k(23)*u(13);
end

