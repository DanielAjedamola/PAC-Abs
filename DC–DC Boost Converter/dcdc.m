%
% dcdc.m
%
% created on: 09.29.2015
%     author: Daniel A.
%
% see readme file for more information on the dcdc example
%
% you need to run ./dcdc binary first 
%
% so that the file: dcdc_controller.bdd is created
%

function dcdc
clear set
close all

%% simulation
% reference initial state
x0 = [1.17 5.6]; %[1.35 5.7];

% number of trajectories to be simulated
N = 1;

% load the symbolic set containing the controller
controller = SymbolicSet('dcdc_controller.bdd');

% simulation horizon
T = 20000;
dt = 0.5;

% % store trajectories in a cell array
% trajectories = cell(N,1);
% 
% for k = 1:N
%     % random initial condition close to x0
%     perturb = 0.01*randn(size(x0)); % small Gaussian noise
%     x_init = x0 + perturb;
% 
%     y = x_init;
%     v = [];
% 
%     for t = 1:T
%         u = controller.getInputs(y(end,:));
%         v = [v; u(1,:)];
%         [~, x] = ode45(@unicycle_ode, [0 dt], y(end,:), ...
%             odeset('abstol',1e-4,'reltol',1e-4), u(1,:));
%         y = [y; x(end,:)];
%     end
% 
%     % store trajectory
%     trajectories{k} = y;
% end

% store trajectories in a cell array
trajectories = cell(N,1);
t_app = 5;

for k = 1:N
    % random initial condition close to x0
    perturb = 0.01*randn(size(x0)); % small Gaussian noise
    x_init = x0 + perturb;

    y = x_init;
    v = [];

    for t = 1:T
        u = controller.getInputs(y(end,:));
        v = [v; u(1,:)];
        [~, x] = ode45(@unicycle_ode, [0 dt], y(end,:), ...
            odeset('abstol',1e-4,'reltol',1e-4), u(1,:));
        y = [y; x(end,:)];
    end

    % store trajectory
    trajectories{k} = y;
end

% eps = 2*0.0005; % approximate value of the ASF, obtained from the .ipynb file
% state_grid = SymbolicSet('dcdc_state_grid.bdd'); % SCOTS grid object
% abs_states = state_grid.points; % all abstract states
% 
% for k = 1:N
%     perturb = 0.01*randn(size(x0));
%     x_init = x0 + perturb;
% 
%     y = x_init;
%     v = [];
% 
%     for t = 1:T
%         current_y = y(end,:);
%         idx = vecnorm(abs_states - current_y, 2, 2) <= eps;
%         u_all = [];
%         for i = find(idx)'
%             u_tmp = controller.getInputs(abs_states(i,:));
%             u_all = [u_all; u_tmp(1,:)];
%         end
% 
%         if ~isempty(u_all)
% %             % Pick a random row (control input) from u_all
% %             idx_rand = randi(size(u_all,1));
% %             u_sel = u_all(idx_rand,:);
% %         else
%             u_sel = controller.getInputs(current_y);
%             u_sel = u_sel(1,:);
%         end
% 
%         v = [v; u_sel];
%         [~, x] = ode45(@unicycle_ode, [0 dt], current_y, ...
%             odeset('abstol',1e-4,'reltol',1e-4), u_sel);
%         y = [y; x(end,:)];
%     end
% 
%     trajectories{k} = y;
% end

%% plot the dcdc safe set
% colors
colors = get(groot,'DefaultAxesColorOrder');

% plot the domain of the controller
x = controller.points;
plot(x(:,1),x(:,2),'.','color',[0.8 0.8 0.8])
hold on

% plot initial states and trajectories
for k = 1:N
    y = trajectories{k};
    plot(y(:,1), y(:,2), '.-', 'color', colors(mod(k-1, size(colors,1))+1,:), 'markersize', 0.5)
    plot(y(1,1), y(1,2), '.', 'color', colors(5,:), 'markersize', 15)
end

% plot safe set
v = [1.15 5.45; 1.55 5.45; 1.15 5.85; 1.55 5.85 ];
patch('vertices', v, 'faces', [1 2 4 3], ...
      'facecolor', 'none', 'edgec', colors(2,:), 'linew', 1)
hold on

box on
grid on
axis([1.1 1.6 5.4 5.9])

end

function dxdt = unicycle_ode(t,x,u)
  r0 = 1.0;
  vs = 1.0;
  rl = 0.05;
  rc = rl / 10;
  xl = 3.0;
  xc = 70.0;

  if (u == 1)
    A = [ -rl / xl, 0 ;
           0, (-1 / xc) * (1 / (r0 + rc)) ] ;
  else
    A = [ (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))), ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
           5 * (r0 / (r0 + rc)) * (1 / xc), (-1 / xc) * (1 / (r0 + rc)) ];
  end
  b = [(vs / xl) ; 0];

  dxdt = A*x + b;
end
