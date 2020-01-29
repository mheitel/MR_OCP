%CSTR_ORIG_MS solve OCP of CSTR via multiple shooting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               optimal control of a CSTR             %%%%
%%% solve optimal control problem with multiple shooting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

clear variables
close all
clc

import casadi.*

saveSolution = 0;
showFigures = 1;
savePlots = 0; % ensure that showFigures = 1, if savePlots = 1!

% % latex options
% set(groot,'DefaultLineLineWidth',2,'DefaultStairLineWidth',2,'DefaultAxesFontSize',16,...
%     'DefaultTextFontSize',16,'DefaultLegendInterpreter','latex')
set(groot,'DefaultLineLineWidth',1.2,'DefaultStairLineWidth',1.2,...
    'DefaultLegendInterpreter','latex','DefaultTextInterpreter','latex')
% set(groot,'DefaultLegendInterpreter','latex','DefaultTextInterpreter','latex');

% % reset options
% reset(groot)

% define some variables
tf = 500;
nx = 4;
nz = 5;
nu = 2;
x_min = 0;
x_max = 1;
q_min = 0;
q_max = 1.5e-3;
qA_max = 1e-3; %m^3/s
V_max = 0.1; %m^3
A_in = 1; %mol/m^3
A0 = 1e-3; % mol/m^3
% A0 = 0.000901; % ->reduced sol, used for cstr_full_N4000-figures
B0 = 1e-3;
C0 = 0;
% C0 = 0.015986; % -> reduced sol, used for cstr_full_N4000-figures
D0 = 1e-8;
V0 = 1e-2; % m^3
V_tf = 1e-2;
epsilon = 1e-6;
k1 = 100; %reaction rate
k1r = 90;
k2 = epsilon;
k3 = 20;
alpha = 1;

% Declare model variables
A = SX.sym('A');
B = SX.sym('B');
C = SX.sym('C');
D = SX.sym('D');
V = SX.sym('V');
z =[A;B;C;D;V];
qA = SX.sym('qA');
q = SX.sym('q');
u = [qA;q];


% Model equations
zdot = [-k1*A+k1r*B+qA/V*(A_in-A);
    k1*A-(k1r+k2)*B-qA/V*B;
    k2*B-k3*C-qA/V*C;
    k3*C-qA/V*D;
    qA-q];

% Continuous time dynamics
f = Function('f', {z, u}, {zdot});

% Control discretization
N = 140; % number of control intervals
M = 1; % RK4 steps per interval
DT = tf/N/M;
integrator = simpleIRK(f,M,3,'radau','newton');
% integrator = simpleRK(f,M,4);

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];


Z0 = MX.sym('Z0', nz);
w = [w, {Z0}];
lbw = [lbw; A0; B0; C0; D0; V0];
ubw = [ubw; A0; B0; C0; D0; V0];
w0  = [w0;  A0; B0; C0; D0; V0];

% Formulate the NLP
Zk = Z0;
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)], nu);
    w = [w, {Uk}];
    lbw = [lbw; q_min*ones(nu,1)];
    ubw = [ubw; qA_max; q_max];
    w0 = [w0; 0.5*(q_min+qA_max);0.5*(q_min+q_max)];
    
    % Integrate till the end of the interval
    % obj = int(-10*q*B+qA^2+q^2,0,tf)
    Zk_end = integrator(Zk, Uk, DT);
    J = J + DT*(-alpha*Uk(2)*Zk(2)+(Uk(1)^2+Uk(2)^2));
    
    % New NLP variable for state at end of interval
    Zk = MX.sym(['Z_' num2str(k+1)], nz);
    w = [w, {Zk}];
    if k<N-1
        lbw = [lbw; x_min*ones(nx,1); 0];
        ubw = [ubw; x_max*ones(nx,1); V_max];
        w0 = [w0; 0.5*(lbw(end-nz+1:end)+ubw(end-nz+1:end))];
    else
        lbw = [lbw; x_min*ones(nx,1); V_tf];
        ubw = [ubw; x_max*ones(nx,1); V_tf];
        w0  = [w0;  0.5*(x_min+x_max)*ones(nx,1); V_tf];
    end
    
    % Add equality constraint
    g = [g, {Zk_end-Zk}];
    lbg = [lbg; zeros(nz,1)];
    ubg = [ubg; zeros(nz,1)];
end

% Create an NLP solver
%opts = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1)); 
opts = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1)); 
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%% Solve the NLP
tic
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
time_full = toc;
w_opt = full(sol.x);

fprintf('elapsed time in total is \t%10.6f seconds\n',time_full);
fprintf('optimal function value: \t%10.6f\n',full(sol.f));

%% Plot the solution
A_opt = w_opt(1:nz+nu:end);
B_opt = w_opt(2:nz+nu:end);
C_opt = w_opt(3:nz+nu:end);
D_opt = w_opt(4:nz+nu:end);
V_opt = w_opt(5:nz+nu:end);
qA_opt = w_opt(nz+1:nz+nu:end);
q_opt = w_opt(nz+2:nz+nu:end);
tgrid = linspace(0, tf, N+1);

if showFigures
    clf;
    figure(1)
    hold on
    
    plot(tgrid, A_opt, tgrid, B_opt,'--',tgrid, 1e7*C_opt,'-.',tgrid, 1e4*D_opt,':')
    legend({'$A$','$B$','$10^{7}\cdot C$','$10^{4}\cdot D$'})
    xlabel('time $t$')
%     ylabel('states')
%     title(['chemical reactor in interval [0,' num2str(tf) '] sec with ' num2str(N) ' subintervals']);
    if savePlots
        saveas(gcf,'../figures/cstrFullStates','epsc')
    end
    
    figure(2)
    hold on
    stairs(tgrid, [qA_opt; nan],'--');
    stairs(tgrid, [q_opt; nan],'--')
    legend({'$q_{\textrm{A}}$','$q$'},'location','East')
    xlabel('time $t$')
%     ylabel('controls')
%     title(['chemical reactor in interval [0,' num2str(tf) '] sec with ' num2str(N) ' subintervals']);
    if savePlots
        saveas(gcf,'../figures/cstrFullControls','epsc')
    end
    
    figure(3)
    hold on
    plot(tgrid, V_opt);
    legend({'$V$'},'location','NorthWest')
    xlabel('time $t$')
%     ylabel('volume')
%     title(['chemical reactor in interval [0,' num2str(tf) '] sec with ' num2str(N) ' subintervals']);
    if savePlots
        saveas(gcf,'../figures/cstrFullVolume','epsc')
    end
       
end

%% use values for initialization of other problem
if saveSolution
    A_full = A_opt;
    B_full = B_opt;
    C_full = C_opt;
    D_full = D_opt;
    V_full = V_opt;
    u1_full = qA_opt;
    u2_full = q_opt;
    save('comparison_full_red/cstr_full_N140.mat','A_full','B_full','C_full','D_full','V_full','u1_full','u2_full','time_full')
    clear A_full B_full C_full D_full V_full u1_full u2_full
end