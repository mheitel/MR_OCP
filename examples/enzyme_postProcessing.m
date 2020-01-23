function [] = enzyme_postProcessing(sol,param, varargin)
%ENZYM_POSTPROCESSING post processing of enzyme example
%   plots the results of the enzyme OCP and saves them (optional) 
% AUTHOR:   Marcus Heitel
% DATE:     Jan 23rd, 2017

import casadi.*
if isempty(sol)
   fprintf(2,'postProcessing not possible. Variable >sol< is empty\n'); 
   return;
end
w_opt = full(sol.x);

%% Plot the solution
x_opt = w_opt(1:3:end);
y_opt = w_opt(2:3:end);
u_opt = w_opt(3:3:end);
tgrid = linspace(0, param.T, param.N+1);
clf;
hold on
plot(tgrid, x_opt, '--')
plot(tgrid, y_opt, '-')
stairs(tgrid, [u_opt; nan], '-.')
xlabel('t')
legend('x','y','u')

% OPTIONAL: save solution in .mat-file 
if nargin>2
    options = varargin{1};
    if options.savesolution
        % generate variables according to suffix
        allvars = {};
        for i = 1:param.nz % states
            name = [char('x'-1+i) '_' suffix];
            eval([name '=' char('x'-1+i) '_opt;']);
            allvars = {allvars{:} name};
        end
        for i = 1:param.nu+1 % time + controls
            name = [char('t'-1+i) '_' suffix];
            eval([name '=' char('t'-1+i) '_opt;']);
            allvars = {allvars{:} name};
        end
        % save variables in savefile
        if options.append
            save([savefile '.mat'],allvars{:} ,'-append')
        else
            save([savefile '.mat'],allvars{:})
        end
    end
end
end