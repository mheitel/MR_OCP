%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Michealis-Menten-Henri enyzme example %%%%
%%% solve original optimal control problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
clear variables
clc

import casadi.*

% Time horizon
T = 5;
epsilon = 1e-6;
eta = 1;

% Declare model variables
x = MX.sym('x');
y = MX.sym('y');
z = [x; y];
u = MX.sym('u');

% dynamics
zdot = [-x+(x+0.5)*y + u; (x-(x+1.0)*y)/epsilon];
f = Function('f', {z, u}, {zdot});
f = f.expand();

% Control discretization
N = 40; % number of control intervals
M = 1; % RK4 steps per interval
DT = T/N/M;
F = simpleIRK(f,M,2,'radau','newton');
% F = simpleRK(fsimple,M,4);

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];
h = {};

% "Lift" initial conditions
Z0 = MX.sym('Z0', 2);
w = [w, {Z0}];
valy = eta/(1+eta) + 0.5*eta./((eta+1).^4)*epsilon;
lbw = [lbw; 0; 0];
ubw = [ubw; 9; 9];
w0 = [w0; eta; valy];
g = [g; {Z0-[eta;valy]}];
lbg = [lbg; 0;0];
ubg = [ubg; 0; 0];
% box constraints for initial conditions
% lbw = [lbw; eta; valy];
% ubw = [ubw; eta; valy];
% w0 = [w0; eta; valy];

% Formulate the NLP
Zk = Z0;
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = [w, {Uk}];
    lbw = [lbw; 0];
    ubw = [ubw; 9];
    w0 = [w0; 1];

    % Integrate till the end of the interval
    Zk_end = F(Zk, Uk, DT);
    J = J + DT*(-50*Zk(2) + Uk^2);

    % New NLP variable for state at end of interval
    Zk = MX.sym(['Z_' num2str(k+1)], 2);
    w = [w, {Zk}];
    lbw = [lbw; 0; 0];
    ubw = [ubw; 9; 9];
    w0 = [w0; eta; 0.5];

    % Add equality constraint
    g = [g, {Zk_end-Zk}];
    lbg = [lbg; 0; 0];
    ubg = [ubg; 0; 0];            
end
h = [h; {vertcat(w{:})-ubw;-vertcat(w{:})+lbw}];

% Create an NLP solver
% opts = struct();
% opts = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1)); 
opts = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1));
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%disp('**** problem within casadi context defined ****')

%% Solve the NLP
tic
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,'lbg', lbg, 'ubg', ubg);
t_opt = toc;
w_opt = full(sol.x);

%disp('**** work is done! ****')
fprintf('elapsed time in total is \t%10.6f seconds\n',t_opt);

%% Plot the solution
x_opt = w_opt(1:3:end);
y_opt = w_opt(2:3:end);
u_opt = w_opt(3:3:end);
tgrid = linspace(0, T, N+1);
clf;
hold on
plot(tgrid, x_opt, '--')
plot(tgrid, y_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x','y','u')