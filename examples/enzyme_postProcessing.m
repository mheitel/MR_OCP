function [] = enzyme_postProcessing(sol,param, varargin)
%ENZYM_POSTPROCESSING post processing of enzyme example
%   plots the results of the enzyme OCP and saves them (optional) 
%
% INPUT: 
%   sol        - solution struct (output of PERFORMREDUCTION) 
%   param      - see input of PERFORMREDUCTION
%   varargin   - contains options for saving results
%    -> savesolution - save solution variables in mat file? (type:bool)
%    -> savefile     - filepath for the saved solution
%                      (only needed if savesolution = true)
%    -> suffix       - suffix of saved filename
%    -> append       - append variables in file? (type:bool)
%
% See also PERFORMREDUCTION
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

import casadi.*
if isempty(sol)
   fprintf(2,'postProcessing not possible. Variable >sol< is empty\n'); 
   return;
end
w_opt = full(sol.x);
% latex options
set(groot,'DefaultLegendInterpreter','latex','DefaultTextInterpreter','latex')

%% Plot the solution
x_opt = w_opt(1:3:end);
y_opt = w_opt(2:3:end);
u1_opt = w_opt(3:3:end);

tgrid = linspace(0, param.T, param.N+1);
clf;
hold on

if nargin>2
    options = varargin{1};
    if isfield(options,'plotForTeX')
        if options.plotForTeX
            reset(groot)
            set(groot,'DefaultLineLineWidth',1.2,...
                'DefaultStairLineWidth',1.2,...
                'DefaultLegendInterpreter','latex',...
                'DefaultTextInterpreter','latex')
        end
    end
end

plot(tgrid, x_opt, '--')
plot(tgrid, y_opt, '-')
stairs(tgrid, [u1_opt; nan], '-.')
xlabel('time $t$')
legend({'$z_s$','$z_f$','$u$'},'location','NorthWest')

% OPTIONAL: save solution in .mat-file 
if nargin>2
    options = varargin{1};
    if options.savesolution
        % generate variables according to suffix
        allvars = {};
        for i = 1:param.nx + param.ny % states
            name = [char('x'-1+i) '_' options.suffix];
            eval([name '=' char('x'-1+i) '_opt;']);
            allvars = {allvars{:} name};
        end
        for i = 1:param.nu % controls
            name = [char('u') num2str(i) '_' options.suffix];
            eval([name '=' char('u') num2str(i) '_opt;']);
            allvars = {allvars{:} name};
        end
        % save time
        name = ['time' '_' options.suffix];
        eval([name '= options.time;'])
        allvars = {allvars{:} name};
        % save variables in savefile
        if options.append
            save([options.savefile '.mat'],allvars{:} ,'-append')
        else
            save([options.savefile '.mat'],allvars{:})
        end
    end
end
end
