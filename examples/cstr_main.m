%CSTR_MAIN solve OCP of CSTR via manifold-based model reduction methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       optimal control of a CSTR      %%%%
%%% solve reduced optimal control problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

clear variables
close all
clc

import casadi.*

% add path of performReduction method
addpath('../src')

%% define parameters for optimization
% number of optimization variables
param.nx = 3;
param.ny = 2;
param.nu = 2;
n_var = param.nx + param.ny + param.nu;
% time scale parameter
% order of ZDP
param.zdp_order = 2;
param.zdp_factor = 1e-3;
% final time
param.T = 500;
% Control discretization
param.N = 140; % number of control intervals
param.M = 1; % RK4 steps per interval

%% initial values and bounds of optimization variables
% (x_1,y_1,u_1,x_2,...,u_n,x_n+1,y_n+1)
q_min = 0;
q_max = 1.5e-3;
qA_max = 1e-3; %m^3/s
V_max = 0.1; %m^3
A_in = 1; %mol/m^3
A0 = 1e-3; % mol/m^3
% A0 = 0.000901;
B0 = 1e-3;
C0 = 0;
D0 = 1e-8;
V0 = 1e-2; % m^3
V_T = 1e-2;
param.w0 = 0.5*ones(n_var*(param.N+1)-param.nu,1);
% volume
param.w0(3:n_var:end) = 0.5*V_max*ones(param.N+1,1);
% param.w0(3:n_var:end) = 0.5*ones(param.N+1,1);
param.w0(1:param.nx) = [B0;D0;V0];
param.w0(param.nx+1:param.nx+param.ny) = [A0;C0];
param.w0(end-param.ny) = V_T;

param.lbw = zeros(n_var*(param.N+1)-param.nu,1);
param.ubw = ones(n_var*(param.N+1)-param.nu,1);
% volume
param.ubw(3:n_var:end) = V_max*ones(param.N+1,1);
% controls
for i=0:param.N-1
   param.lbw(param.nx+param.ny+1+i*n_var:(i+1)*n_var) = [q_min;q_min]; 
   param.ubw(param.nx+param.ny+1+i*n_var:(i+1)*n_var) = [qA_max;q_max]; 
   param.w0(param.nx+param.ny+1+i*n_var:(i+1)*n_var) = 0.5*[q_min;q_min]+0.5*[qA_max;q_max]; 
end

param.lbw(1:param.nx) = [B0;D0;V0];
param.ubw(1:param.nx) = [B0;D0;V0];
param.lbw(end-param.ny) = V_T;
param.ubw(end-param.ny) = V_T;

param.lbw(param.nx+2) = C0;
param.ubw(param.nx+2) = C0;

%% define dynamics (ode +  objective function)
% obective function/integrand
k1 = 100; %reaction rate
k1r = 90;
k2 = 1e-6;%param.epsilon;
k3 = 20;
alpha = 1;
dynamics.L = @(z,u) -alpha*u(2)*z(1)+(u(1)^2+u(2)^2);
% dynamics - ode
dynamics.x = @(z,u) [k1*z(4)-(k1r+k2)*z(1)-u(1)/z(3)*z(1);
    k3*z(5)-u(1)/z(3)*z(2);
    u(1)-u(2)];
dynamics.y = @(z,u) [-k1*z(4)+k1r*z(1)+u(1)/z(3)*(A_in-z(4));
    k2*z(1)-k3*z(5)-u(1)/z(3)*z(5)];

%% define options (not neccessary)
options.savesolution = false;
options.saveplots = false;
options.savefile = ['./path/to/folder/cstr_lift_zdp_N' num2str(param.N)];
options.suffix = 'lift_zdp';
options.append = false;
%options.ipopt_options = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1));
options.ipopt_options = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1));
options.plotForTeX = true;

%% perform optimization
% call template for model reduction
[sol, time] = performReduction(dynamics, param, options);
% post Processing -> Plot, save files
options.time = time;
cstr_postProcessing(sol,param, dynamics, options)
