%ENZYME_ORIG_CO solve enzyme OCP via collocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Michealis-Menten-Henri enyzme example     %%%%
%%% solve optimal control problem with collocation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

clear variables
close all
%clc

import casadi.*

saveSolution = 0;

% Degree of interpolating polynomial
d = 2;

% Get collocation points
tau_root = [0 collocation_points(d, 'radau')];

[B,C,D] = collocation(tau_root);

% Time horizon
T = 5;
epsilon = 1e-6;
eta = 1;

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x = [x1; x2];
u = SX.sym('u');

% Model equations
xdot = [-x1+(x1+0.5)*x2 + u; (x1-(x1+1.0)*x2)/epsilon];

% Objective term
L = -50*x2 + u^2;

% Continuous time dynamics
f = Function('f', {x, u}, {xdot, L});

% Control discretization
N = 40; % number of control intervals
h = T/N;
dimx = 2;
dimu = 1;

% Start with an empty NLP
w=cell(1,N*(2+d)+1);
w_ind = 1;
w0 = zeros(dimx*(N*(d+1)+1)+dimu*N,1);
lbw = zeros(dimx*(N*(d+1)+1)+dimu*N,1);
ubw = zeros(dimx*(N*(d+1)+1)+dimu*N,1);
wb_ind = 1;
J = 0;
g=cell(1,N*(d+1));
g_ind = 1;
lbg = zeros(dimx*N*(d+1),1);
ubg = zeros(dimx*N*(d+1),1);
gb_ind = 1;

% "Lift" initial conditions
X0 = SX.sym('X0', 2);
w(w_ind)= {X0};w_ind = w_ind+1;
lbw(wb_ind:wb_ind+dimx-1) = [eta; eta/(1.0+eta)];
ubw(wb_ind:wb_ind+dimx-1) = [eta; eta/(1.0+eta)];
w0(wb_ind:wb_ind+dimx-1) = [eta; eta/(1.0+eta)];
wb_ind = wb_ind + dimx;

%% Formulate the NLP
Xk = X0;
for k=0:N-1
    % New NLP variable for the control
    Uk = SX.sym(['U_' num2str(k)]);
    w(w_ind) = {Uk};w_ind=w_ind+1;
    lbw(wb_ind:wb_ind+dimu-1) = 0;
    ubw(wb_ind:wb_ind+dimu-1) = 9;
    w0(wb_ind:wb_ind+dimu-1) = 1;
    wb_ind = wb_ind + dimu;
    
    % State at collocation points
    Xkj = cell(d,1);
    for j=1:d
        Xkj{j} = SX.sym(['X_' num2str(k) '_' num2str(j)], 2);
        w(w_ind) = Xkj(j);
        w_ind = w_ind+1;
        lbw(wb_ind:wb_ind+dimx-1) = [0;0];
        ubw(wb_ind:wb_ind+dimx-1) = [9;9];
        w0(wb_ind:wb_ind+dimx-1) = [4;1];
        wb_ind = wb_ind + dimx;
    end
    
    % Loop over collocation points
    Xk_end = D(1)*Xk;
    for j=1:d
        % Expression for the state derivative at the collocation point
        xp = C(1,j+1)*Xk;
        for r=1:d
            xp = xp + C(r+1,j+1)*Xkj{r};
        end
        
        % Append collocation equations
        [fj, qj] = f(Xkj{j},Uk);
        g(g_ind) = {h*fj-xp};
        g_ind = g_ind + 1;
        lbg(gb_ind:gb_ind+dimx-1) = [0; 0];
        ubg(gb_ind:gb_ind+dimx-1) = [0; 0];
        gb_ind = gb_ind  + dimx;
        
        % Add contribution to the end state
        Xk_end = Xk_end + D(j+1)*Xkj{j};
        
        % Add contribution to quadrature function
        J = J + B(j+1)*qj*h;
    end
    
    % New NLP variable for state at end of interval
    Xk = SX.sym(['X_' num2str(k+1)], 2);
    w(w_ind) = {Xk};
    w_ind = w_ind +1 ;
    lbw(wb_ind:wb_ind+dimx-1) = [0;0];
    ubw(wb_ind:wb_ind+dimx-1) = [9;9];
    w0(wb_ind:wb_ind+dimx-1) = [1;1];
    wb_ind = wb_ind + dimx;
    
    % Add equality constraint
    g(g_ind)  = {Xk_end-Xk};
    g_ind = g_ind + 1;
    lbg(gb_ind:gb_ind+dimx-1) = [0; 0];
    ubg(gb_ind:gb_ind+dimx-1) = [0; 0];
    gb_ind = gb_ind  + dimx;
end

%% Create an NLP solver
% opts = struct();
%opts = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1)); 
opts = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1)); 
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%% Solve the NLP
tic
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
time_co_full = toc;
w_opt = full(sol.x);
fprintf('elapsed time in total is \t%10.6f seconds\n',time_co_full);

%% Plot the solution
x1_opt = w_opt(1:3+2*d:end);
x2_opt = w_opt(2:3+2*d:end);
u_opt = w_opt(3:3+2*d:end);
tgrid = linspace(0, T, N+1);
figure
hold on
plot(tgrid, x1_opt, '--')
plot(tgrid, x2_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x1','x2','u')

%% save solution
if saveSolution
    x_co_full = x1_opt;
    y_co_full = x2_opt;
    u1_co_full = u_opt;  
    save('comparison_full_red/enzym_co_full_N40.mat','x_co_full', 'y_co_full', 'u1_co_full' ,'time_co_full'); % '-append'
    clear x_co_full y_co_full u1_co_full   
end