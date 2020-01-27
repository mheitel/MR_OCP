function [] = cstr_postProcessing(sol,param, dynamics, varargin)

%               -> savesolution - save solution variables in mat file?
%                                 (type:bool)
%               -> savefile     - filepath for the saved solution
%                                 (only needed if savesolution = true)
%               -> suffix       -
%               -> append       -
if isempty(sol)
    fprintf(2,'postProcessing not possible. Variable >sol< is empty\n');
    return;
end
import casadi.*

w_opt = full(sol.x);
nz = param.nx + param.ny;
N = param.N;

%% Plot the solution
q1_opt = w_opt(nz+1:nz+param.nu:end);
q2_opt = w_opt(nz+2:nz+param.nu:end);

%% TEMP TEST
t_fine = linspace(0,param.T,N+1);
AAA = w_opt(4:nz+param.nu:end);
BBB = w_opt(1:nz+param.nu:end);
CCC = w_opt(5:nz+param.nu:end);
DDD = w_opt(2:nz+param.nu:end);
VVV = w_opt(3:nz+param.nu:end);
set(groot,'DefaultLegendInterpreter','latex','DefaultTextInterpreter','latex');
plot(t_fine, AAA, t_fine, BBB,'--',t_fine, 1e7*CCC ,'-.',t_fine, 1e4*DDD,':',t_fine,VVV)
legend({'$A$','$B$','$10^{7}\cdot C$','$10^{4}\cdot D$','$V$'},'location','East')
xlabel('time $t$')
%%

% fprintf('D opt: %.4e\n',w_opt(2:nz+param.nu:end))
tgrid = linspace(0, param.T, N+1);


%% integrating forward with optimal control
% figure
% n_int = 10*N;
n_int = N;
dt = param.T/n_int;
t_fine = linspace(0,param.T,n_int+1);
Z = zeros(nz,n_int);
Z(:,1) = w_opt(1:nz);
% integrator
z = MX.sym('z',nz);
u = MX.sym('u',param.nu);
xdot = dynamics.x(z,u);
ydot = dynamics.y(z,u);
obj = Function('obj',{z,u},{dynamics.L(z,u)});
f_tilde = Function('f_tilde',{z,u},{[xdot;ydot]});
F = simpleIRK(f_tilde,param.M,3,'radau','newton');
obj_approx = 0;
for i=1:n_int
    Z(:,i+1) = full(F(Z(:,i),[q1_opt(ceil(N*i/n_int));q2_opt(ceil(N*i/n_int))],dt));
    obj_approx = obj_approx + dt*obj(Z(:,i),[q1_opt(ceil(N*i/n_int));q2_opt(ceil(N*i/n_int))]);
end
fprintf('approx optimal function value: \t\t%10.16f\n',full(obj_approx));

%% fiugures
if nargin > 3
    options = varargin{1};
else
    options.saveplots = false;
end
alreadyPlotted = false;

if nargin>2
    options = varargin{1};
    if isfield(options,'plotForTeX')
        if options.plotForTeX
            set(groot,'DefaultLineLineWidth',1.2,...
                'DefaultStairLineWidth',1.2,...
                'DefaultLegendInterpreter','latex',...
                'DefaultTextInterpreter','latex')
            
            %plot used in LaTeX
            figure
            hold on
            plot(t_fine, Z(4,:), t_fine, Z(1,:),'--',t_fine, 1e7*Z(5,:),'-.',t_fine, 1e4*Z(2,:),':')
            legend({'$A$','$B$','$10^{7}\cdot C$','$10^{4}\cdot D$'},'location','East')
            xlabel('time $t$')
            % title(['chemical reactor in interval [0,' num2str(param.T) '] sec with ' num2str(N) ' subintervals']);
            if options.saveplots
                saveas(gcf,'../figures/cstrZDPStates','epsc')
            end
            
            figure
            hold on
%             stairs(tgrid, [q1_opt; nan],'-');
%             stairs(tgrid, [q2_opt; nan],'--')
            plot(tgrid, [q1_opt; nan],'-');
            plot(tgrid, [q2_opt; nan],'--')
            legend({'$q_{\textrm{A}}$','$q$'},'location','East')
            xlabel('time $t$')
            % title(['chemical reactor in interval [0,' num2str(param.T) '] sec with ' num2str(N) ' subintervals']);
            if options.saveplots
                saveas(gcf,'../figures/cstrZDPControls','epsc')
            end
            
            figure
            hold on
            plot(t_fine, Z(3,:));
            legend({'$V$'},'location','NorthWest')
            xlabel('time $t$')
            % title(['chemical reactor in interval [0,' num2str(param.T) '] sec with ' num2str(N) ' subintervals']);
            if options.saveplots
                saveas(gcf,'../figures/cstrZDPVolume','epsc')
            end
            
            alreadyPlotted = true;
        end
    end
end
if (~alreadyPlotted)
    set(groot,'DefaultLegendInterpreter','latex','DefaultTextInterpreter','latex');
    
    % standard plots
    figure
    hold on
    plot(t_fine, Z(4,:), t_fine, Z(1,:),'--',t_fine, 1e7*Z(5,:),'-.',t_fine, 1e4*Z(2,:),':')
    legend({'$A$','$B$','$10^{7}\cdot C$','$10^{4}\cdot D$'},'location','East')
    xlabel('time $t$')
    ylabel('states')
    % title(['chemical reactor in interval [0,' num2str(param.T) '] sec with ' num2str(N) ' subintervals']);
    if options.saveplots
        saveas(gcf,'../figures/cstrZDPStates','epsc')
    end
    
    figure
    hold on
    stairs(tgrid, [q1_opt; nan],'-');
    stairs(tgrid, [q2_opt; nan],'--')
    legend({'$q_{\textrm{A}}$','$q$'},'location','East')
    xlabel('time $t$')
%     ylabel('controls')
    % title(['chemical reactor in interval [0,' num2str(param.T) '] sec with ' num2str(N) ' subintervals']);
    if options.saveplots
        saveas(gcf,'../figures/cstrZDPControls','epsc')
    end
    
    figure
    hold on
    plot(t_fine, Z(3,:));
    legend({'$V$'},'location','NorthWest')
    xlabel('time $t$')
%     ylabel('volume')
    % title(['chemical reactor in interval [0,' num2str(param.T) '] sec with ' num2str(N) ' subintervals']);
    if options.saveplots
        saveas(gcf,'../figures/cstrZDPVolume','epsc')
    end
end


% OPTIONAL: save solution in .mat-file
if nargin>2
    if options.savesolution
        % generate variables according to suffix
        allvars = {};
        A_opt = Z(4,:)';
        B_opt = Z(1,:)';
        C_opt = Z(5,:)';
        D_opt = Z(2,:)';
        V_opt = Z(3,:)';
        for i = 1:nz-1 % states
            name = [char('A'-1+i) '_' options.suffix];
            eval([name '=' char('A'-1+i) '_opt;']);
            allvars = {allvars{:} name};
        end
        name = ['V_' options.suffix];
        eval([name '= V_opt;']);
        allvars = {allvars{:} ['V_' options.suffix]};
        for i = 1:param.nu % time + controls
            name = [char('u') num2str(i) '_' options.suffix];
            eval([name '=' char('q') num2str(i) '_opt;']);
            allvars = {allvars{:} name};
        end
        % save computing time
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
    

% reset options
reset(groot)
end
