%ENZYME_ORDER_RED solve enzyme OCP via order reduced model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          Michealis-Menten-Henri enyzme example      %%%%
%%% solve order reduced problem with ZDP approximation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

clear variables
close all
%clc

import casadi.*

% Time horizon
T = 5;
epsilon = 1e-6;
eta = 1;

% parameters for solving bvp in order to calculate slow manifold
d2 = 2;% order of interpolating polynomial
N2 = 1;% number of subintervals of bvp

% Declare model variables
nx = 1;
ny = 1;
nz = nx+ny;
nu = 1;
x = MX.sym('x',nx);
u = MX.sym('u',nu);
[sim, offset] = enzyme_sim_zdp(epsilon);
tmp = sim(x,u);
y = tmp(offset+1:offset+ny);

% Model equations
xdot = -x+(x+0.5)*y + u;

% Continuous time dynamics
f = Function('f', {x, u}, {xdot});

% Control discretization
N = 40; % number of control intervals
M = 1; % RK4 steps per interval
DT = T/N/M;
F = simpleRK(f,M,4);

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
X0 = MX.sym('X0', nx);
w = [w, {X0}];
lbw = [lbw; eta];
ubw = [ubw; eta];
w0 = [w0; eta];

% Formulate the NLP
Xk = X0;
Uk = MX.sym('U_0', nu);
w = [w, {Uk}];
lbw = [lbw; 0];
ubw = [ubw; 9];
w0 = [w0; 1];

for k=0:N-1
    if k>0
        % New NLP variable for the control
        Uk = MX.sym(['U_' num2str(k)], nu);
        w = [w, {Uk}];
        lbw = [lbw; 0];
        ubw = [ubw; 9];
        w0 = [w0; 1];
    end
    
    % get fast variables
        tmp = sim(Xk,Uk);
        Yk = tmp(offset+1:offset+ny);
    
    % Integrate till the end of the interval
    Xk_end = F(Xk, Uk, DT);
    J = J + DT*(-50*Yk + Uk^2);
    
    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], nx);
    w = [w, {Xk}];
    lbw = [lbw; 0];
    ubw = [ubw; 9];
    w0 = [w0; eta];
    
    % Add equality constraint
    g = [g, {Xk_end-Xk}];
    lbg = [lbg; 0];
    ubg = [ubg; 0];
end

% Create an NLP solver
%opts = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1)); 
opts = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1)); 
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%% Solve the NLP
tic
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
t = toc;
w_opt = full(sol.x);
fprintf('elapsed time in total is \t%10.6f seconds\n',t);

%% Plot the solution
x_opt = w_opt(1:nx+nu:end);
u_opt = w_opt(nx+1:nx+nu:end);
y_opt = zeros(length(x_opt),1);
for i=1:length(y_opt)-1
    tmp = full(sim(x_opt(i),u_opt(i)));
    y_opt(i) = tmp(offset+1:offset+ny);
end
tmp = full(sim(x_opt(end),u_opt(end)));
y_opt(end) = tmp(offset+1:offset+ny);
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x_opt, '--')
plot(tgrid, y_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
% legend('x','y','u')