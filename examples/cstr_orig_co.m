%CSTR_ORIG_CO solve OCP of CSTR via collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            optimal control of a CSTR          %%%%
%%% solve optimal control problem with collocation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

set(groot,'DefaultLegendInterpreter','latex', ...
    'DefaultTextInterpreter','latex');

% % reset options
% reset(groot)

% Degree of interpolating polynomial
d = 2;

% Get collocation points
tau_root = [0 collocation_points(d, 'radau')];

[B,C,D] = collocation(tau_root);

% define some variables
tf = 500;
nx = 4;
nz = 5;
nu = 2;
z_min = 0;
z_max = 1;
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
zA = SX.sym('zA');
zB = SX.sym('zB');
zC = SX.sym('zC');
zD = SX.sym('zD');
zV = SX.sym('zV');
z =[zA;zB;zC;zD;zV];
qA = SX.sym('qA');
q = SX.sym('q');
u = [qA;q];


% Model equations
zdot = [-k1*zA+k1r*zB+qA/zV*(A_in-zA);
    k1*zA-(k1r+k2)*zB-qA/zV*zB;
    k2*zB-k3*zC-qA/zV*zC;
    k3*zC-qA/zV*zD;
    qA-q];

% Objective term
L = -alpha*q*zB + qA.^2 + q.^2;

% Continuous time dynamics
f = Function('f', {z, u}, {zdot, L});

% Control discretization
N = 140; % number of control intervals
h = tf/N;
dimz = 5;
dimu = 2;

% Start with an empty NLP
w=cell(1,N*(2+d)+1);
w_ind = 1;
w0 = zeros(dimz*(N*(d+1)+1)+dimu*N,1);
lbw = zeros(dimz*(N*(d+1)+1)+dimu*N,1);
ubw = zeros(dimz*(N*(d+1)+1)+dimu*N,1);
wb_ind = 1;
J = 0;
g=cell(1,N*(d+1));
g_ind = 1;
lbg = zeros(dimz*N*(d+1),1);
ubg = zeros(dimz*N*(d+1),1);
gb_ind = 1;

% "Lift" initial conditions
Z0 = SX.sym('Z0', dimz);
w(w_ind)= {Z0};w_ind = w_ind+1;
lbw(wb_ind:wb_ind+dimz-1) = [A0;B0;C0;D0;V0];
ubw(wb_ind:wb_ind+dimz-1) = [A0;B0;C0;D0;V0];
w0(wb_ind:wb_ind+dimz-1) = [A0;B0;C0;D0;V0];
wb_ind = wb_ind + dimz;

%% Formulate the NLP
Zk = Z0;
for k=0:N-1
    % New NLP variable for the control
    Uk = SX.sym(['U_' num2str(k)],dimu);
    w(w_ind) = {Uk};w_ind=w_ind+1;
    lbw(wb_ind:wb_ind+dimu-1) = q_min*ones(dimu,1);
    ubw(wb_ind:wb_ind+dimu-1) = [qA_max; q_max];
    w0(wb_ind:wb_ind+dimu-1) = [0.5*(q_min+qA_max);0.5*(q_min+q_max)];
    wb_ind = wb_ind + dimu;
    
    % State at collocation points
    Zkj = cell(d,1);
    for j=1:d
        Zkj{j} = SX.sym(['X_' num2str(k) '_' num2str(j)], dimz);
        w(w_ind) = Zkj(j);
        w_ind = w_ind+1;
        lbw(wb_ind:wb_ind+dimz-1) = -inf*ones(dimz,1);
        ubw(wb_ind:wb_ind+dimz-1) = inf*ones(dimz,1);
        w0(wb_ind:wb_ind+dimz-1) = [A0;B0;C0;D0;V0];
        wb_ind = wb_ind + dimz;
    end
    
    % Loop over collocation points
    Zk_end = D(1)*Zk;
    for j=1:d
        % Expression for the state derivative at the collocation point
        zp = C(1,j+1)*Zk;
        for r=1:d
            zp = zp + C(r+1,j+1)*Zkj{r};
        end
        
        % Append collocation equations
        [fj, qj] = f(Zkj{j},Uk);
        g(g_ind) = {h*fj-zp};
        g_ind = g_ind + 1;
        lbg(gb_ind:gb_ind+dimz-1) = zeros(dimz,1);
        ubg(gb_ind:gb_ind+dimz-1) = zeros(dimz,1);
        gb_ind = gb_ind  + dimz;
        
        % Add contribution to the end state
        Zk_end = Zk_end + D(j+1)*Zkj{j};
        
        % Add contribution to quadrature function
        J = J + B(j+1)*qj*h;
    end
    
    % New NLP variable for state at end of interval
    Zk = SX.sym(['Z_' num2str(k+1)], dimz);
    w(w_ind) = {Zk};
    w_ind = w_ind +1 ;
    if k<N-1
        lbw(wb_ind:wb_ind+dimz-1) = [z_min*ones(dimz-1,1);0];
        ubw(wb_ind:wb_ind+dimz-1) = [z_max*ones(dimz-1,1);V_max];
        w0(wb_ind:wb_ind+dimz-1) = 0.5*(lbw(end-nz+1:end)+ubw(end-nz+1:end));
    else
        lbw(wb_ind:wb_ind+dimz-1) = [z_min*ones(dimz-1,1);V_tf];
        ubw(wb_ind:wb_ind+dimz-1) = [z_max*ones(dimz-1,1);V_tf];
        w0(wb_ind:wb_ind+dimz-1) = [0.5*(lbw(end-nz+1:end-1)+ubw(end-nz+1:end-1));V_tf];
    end
    wb_ind = wb_ind + dimz;
    
    % Add equality constraint
    g(g_ind)  = {Zk_end-Zk};
    g_ind = g_ind + 1;
    lbg(gb_ind:gb_ind+dimz-1) = zeros(dimz,1);
    ubg(gb_ind:gb_ind+dimz-1) = zeros(dimz,1);
    gb_ind = gb_ind  + dimz;
end

%% Create an NLP solver
%opts = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1));
opts = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1));
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%% Solve the NLP and measure time
tic
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
t = toc;

w_opt = full(sol.x);

fprintf('elapsed time in total is \t%10.6f seconds\n',t);

%% Plot the solution
var_set = dimz*(d+1) + dimu;
A_opt = w_opt(1:var_set:end);
B_opt = w_opt(2:var_set:end);
C_opt = w_opt(3:var_set:end);
D_opt = w_opt(4:var_set:end);
V_opt = w_opt(5:var_set:end);
qA_opt = w_opt(dimz+1:var_set:end);
q_opt = w_opt(dimz+2:var_set:end);
tgrid = linspace(0, tf, N+1);

if showFigures
    clf;
    figure(1)
    hold on
    
    plot(tgrid, A_opt, tgrid, B_opt,tgrid, 1e7*C_opt,tgrid, 1e4*D_opt)
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
    legend({'$q_{\textrm{A}}$','$q$'})
    xlabel('time $t$')
%     ylabel('controls')
%     title(['chemical reactor in interval [0,' num2str(tf) '] sec with ' num2str(N) ' subintervals']);
    if savePlots
        saveas(gcf,'../figures/cstrFullControls','epsc')
    end
    
    figure(3)
    hold on
    plot(tgrid, V_opt);
    legend({'$V$'})
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
    save('comparison_full_red/cstr_full_N40.mat','A_full','B_full','C_full','D_full','V_full','u1_full','u2_full','time_full')
    clear A_full B_full C_full D_full V_full u1_full u2_full
end
