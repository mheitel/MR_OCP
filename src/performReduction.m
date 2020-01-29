function [sol, time] = performReduction(dynamics,param, varargin)
%PERFORMREDUCTION performs manifold based model reduction of OCP
%   sol = performReduction(dynamics,param, options) attempts to solve the
%   following optimal control problem with singularly perturbed dynamics
%           min int(L,0,T)dt
%           subject to: dx/dt       = f_s(x,y,u)
%                       eps*dy/dt   = f_f(x,y,u)
%                       lbw <= w=(x,y,u) <= ubw
%
% INPUT:    dynamics    - struct for dynamics
%               -> L            - function handle for integrand L
%               -> x            - function handle for ODE in x, i.e. f_s
%               -> y            - function handle for ODE in y, i.e. f_f
%           param       - struct which contains parameters for the OCP
%               -> nx           - number of slow variables x (length(x))
%               -> ny           - number of fast variables y (length(y))
%               -> nu           - number of controls u (length(u))
%               -> zdp_order    - order of time derivative for the Zero
%                                 Derivative Principle
%               -> zdp_factor   - scaling factor for root finding problem
%                                 ZDP (i.e factor*g = 0 instead of g=0)
%                                 e.g. the time scale parameter epsilon
%                                 in case of singular perturbed ODEs.
%               -> T            - end time T of OCP
%               -> N            - number of control intervals (we use
%                                 piecewise constant controls)
%               -> M            - number of integrator steps per interval
%               -> w0           - initial values for optimization variables
%                                 w=(X_0,Y_0,U_0,X_1,...,U_N-1,X_N,Y_N)
%               -> lbw          - lower bound for optimization variables
%               -> ubw          - upper bound for optimization variables
%           options     - optional struct for IPOPT options etc.
%               -> ipopt_options- struct with IPOPT options
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017
%
import casadi.*

% number of variables
nz = param.nx + param.ny;
% step length of discretization
DT = param.T/param.N/param.M;

% Declare model variables
x = SX.sym('x',param.nx);
y = SX.sym('y',param.ny);
z = [x; y];
u = SX.sym('u',param.nu);

% dynamics
xdot = dynamics.x(z,u);
ydot = dynamics.y(z,u);
zdot = [xdot;ydot];
% objective function
obj = Function('obj',{z,u},{dynamics.L(z,u)});

% Continuous time dynamics
f = Function('f', {x,[y;u]}, {xdot});
f = f.expand();
f_tilde = Function('f_tilde',{z,u},{zdot});

% use ZDP for the approximative calculation of the SIM
% switch order (=k-th time-derivative of slow components) of ZDP
if param.zdp_order == 1
    sim_const = Function('sim_const',{z,u}, {param.zdp_factor*ydot});
    sim_const = sim_const.expand();
elseif param.zdp_order == 2
    % ffast = Function('ffast',{y,[x;u]},{ydot});
    jac_ffast = Function('jac_ffast',{z,u},{jacobian(ydot,z)});
    sim_const = Function('sim_const',{z,u}, {param.zdp_factor^2*jac_ffast(z,u)*f_tilde(z,u)});
    sim_const = sim_const.expand();
elseif param.zdp_order == -1
    % case of local BVP approach for SIM instead of ZDP
    jac_f_tilde = Function('jac_f_tilde',{z,u},{jacobian(zdot,z)});
    obj_sim = Function('obj',{z,u}, {0.5*f_tilde(z,u)'*jac_f_tilde(z,u)'*jac_f_tilde(z,u)*f_tilde(z,u)});
    z_feasible = 0.5*param.lbw(nz+param.nu+1:2*nz+param.nu)+0.5*param.ubw(nz+param.nu+1:2*nz+param.nu);
    u_feasible = 0.5*param.lbw(nz+1:nz+param.nu)+0.5*param.ubw(nz+1:nz+param.nu);
    sim_const = Function('sim_const',{z,u},{jacobian(obj_sim(z,u)/obj_sim(z_feasible,u_feasible),y)'});
    sim_const = sim_const.expand();
else
    fprintf(2,'Only ZDP order 1 and 2 are supported -> set zdp_order = 2\n');
    fprintf(2,'ZDP order -1 is the special case of a local BVP approach\n');
    fprintf(2,'Your choice was zdp_order = %d\n', param.zdp_order);
    sol = [];
    return;
end

% integrator = simpleIRK(f,M,2,'radau','newton');
integrator = simpleRK(f,param.M,4); % RK4

% Start with an empty NLP
w={};
J = 0;
g={};
lbg = [];
ubg = [];

% Initial conditions
Z0 = SX.sym('Z0', nz);
w = [w, {Z0}];

% New NLP variable for the control
Uk = SX.sym('U_0',param.nu);
w = [w, {Uk}];

% % Link fast and slow components via manifold
% g = [g,{sim_const(Z0,Uk)}];
% lbg = [lbg; zeros(param.ny,1)];
% ubg = [ubg; zeros(param.ny,1)];

% Formulate the NLP
Zk = Z0;
for k=0:param.N-1
    % Integrate till the end of the interval
    Xk_end = integrator(Zk(1:param.nx),[Zk(param.nx+1:nz); Uk], DT); % only for slow variables
    % piecewise constant approximation of objective integral
    J = J + DT*obj(Zk,Uk);
    
    % New NLP variable for state at the end of the interval
    Zk = SX.sym(['Z_' num2str(k+1)], nz);
    w = [w, {Zk}];
    
    % Add equality constraint(s)
    g = [g, {Xk_end-Zk(1:param.nx)}];
    lbg = [lbg; zeros(param.nx,1)];
    ubg = [ubg; zeros(param.nx,1)];
    
    if k<param.N-1
        % New NLP variable for the control
        Uk = SX.sym(['U_' num2str(k)],param.nu);
        w = [w, {Uk}];
    end
    
    % Link fast and slow components via manifold
    g = [g, {sim_const(Zk,Uk)}];
    lbg = [lbg; zeros(param.ny,1)];
    ubg = [ubg; zeros(param.ny,1)];
end

% Create an NLP solver -> IPOPT
opts = struct();
if nargin > 2 && isfield(varargin{1},'ipopt_options')
    opts = varargin{1}.ipopt_options;
end
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob,opts);

%% Solve the NLP
tic
sol = solver('x0', param.w0, 'lbx', param.lbw, 'ubx', param.ubw,'lbg', lbg, 'ubg', ubg);
time = toc;
fprintf('elapsed time in total is \t%10.6f seconds\n',time);
% fprintf('optimal function value: \t%10.6f\n',full(sol.f));
