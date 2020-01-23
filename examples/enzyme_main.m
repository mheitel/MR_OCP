%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Michealis-Menten-Henri enyzme example %%%%
%%% solve reduced optimal control problem %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

clear variables
close all
clc

% add path of performReduction method
addpath('../src')

%% define parameters for optimization
% number of optimization variables
param.nx = 1;
param.ny = 1;
param.nu = 1;
n_var = param.nx + param.ny + param.nu;
% time scale parameter
epsilon = 1e-6;
% order of ZDP
param.zdp_order = 2;
param.zdp_factor = epsilon;
% final time
param.T = 5;
% Control discretization
param.N = 40; % number of control intervals
param.M = 1; % RK4 steps per interval
% fixed initial value for slow components
x0 = 1;
% initial values of optimization variables
% (x_1,y_1,u_1,x_2,...,u_n,x_n+1,y_n+1)
param.w0 = ones(n_var*(param.N+1)-param.nu,1);
param.w0(1:n_var:end) = x0*ones(param.N+1,1);
param.w0(2:n_var:end) = 0.5*ones(param.N+1,1);
% lower and upper bounds for optimization variables
param.lbw = zeros(n_var*(param.N+1)-param.nu,1);
param.lbw(3:n_var:end) = 0*ones(param.N,1);
param.ubw = 9*ones(n_var*(param.N+1)-param.nu,1);
param.ubw(3:n_var:end) = 9*ones(param.N,1);
param.lbw(1) = x0;
param.ubw(1) = x0;
param.w0(1) = x0;

%% define dynamics (ode +  objective function)
% obective function/integrand
dynamics.L = @(z,u) -50*z(2) + u(1).^2;
% dynamics - ode
dynamics.x = @(z,u) -z(1)+(z(1)+0.5)*z(2) + u(1);
dynamics.y = @(z,u) (z(1)-(z(1)+1.0)*z(2))/epsilon;

%% define options (not neccessary)
options.savesolution = false;
options.savefile = ['./path/to/folder/enzym_N' num2str(param.N)];
options.suffix = 'lift_zdp2';
options.append = false;
%options.ipopt_options = struct('print_time',false,'ipopt',struct('linear_solver','ma27','print_level',1));
options.ipopt_options = struct('print_time',false,'ipopt',struct('linear_solver','mumps','print_level',1));

%% perform optimization
% call template for model reduction
sol = performReduction(dynamics, param, options);
% post Processing -> Plot, save files
enzyme_postProcessing(sol,param, options)
