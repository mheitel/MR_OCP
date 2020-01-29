function [y_sim,offset] = enzyme_sim_zdp(epsilon)
% ENZYME_SIM_ZDP returns function that calculates point of SIM for michaelis
% menten mechanism using the ROOTFINDER function of casadi
% it solves the problem:
%       d^2 g/dt^2 (rpv,y)= 0 -> y 
%   INPUT: 
%       - epsilon is a parameter for the time scale separation
%   OUTPUT
%       - y_sim = h that maps rpv=x(t0), u to x(t0),y(t0) 
%       - offset -index offset for y in output y_sim
%
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

import casadi.*

%% Problem definition
nx = 1;
ny = 1;
nu = 1;
x = MX.sym('x',nx);
y = MX.sym('y',ny);
z = [x;y];
u = MX.sym('u',nu);
zdp_factor = epsilon;

ydot = (x-(x+1.0)*y)/epsilon;
zdot = [-x+(x+0.5)*y + u; ydot];
f_tilde = Function('f_tilde',{z,u},{zdot});
jac_ffast = Function('jac_ffast',{z,u},{jacobian(ydot,z)});
sim_const = Function('sim_const',{y,x,u}, {zdp_factor^2*jac_ffast(z,u)*f_tilde(z,u)});
sim_const = sim_const.expand();

%% create rootfinder function (of newton-type)
rf = rootfinder('rf','newton',sim_const,struct('verbose',false));
y_sim = Function('y_sim',{x,u},{rf(0.5*ones(ny,1),x,u)}); % zones as initial values for lam
offset = 0;
