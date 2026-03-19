%
% engine.m
%
% created on: Nov.05.2015
%     author: Daniel A.
%
% see readme file for more information on the engine example
%
% you need to run ./engine binary first 
%
% so that the file: engine_controller.bdd is created
%

function engine
clear set
close all

%% simulation
% reference initial state
% x0 = [0.4489 0.6483];
% Define intervals where for random initialization
I1 = [0.4488, 0.4490];   % interval for x(1)
I2 = [0.6482, 0.6484];   % interval for x(2)

% number of trajectories to be simulated
N = 10;

% load the symbolic set containing the controller
controller = SymbolicSet('engine_controller.bdd');

% simulation horizon
T = 1000;
dt = 0.1;

% store trajectories in a cell array
trajectories = cell(N,1);

for k = 1:N
%     % random initial condition close to x0
%     perturb = 0.00001*randn(size(x0)); % small Gaussian noise
%     x_init = x0 + perturb;
%     Sample initial state uniformly from I1 × I2
    x_init = [I1(1) + (I1(2) - I1(1)) * rand() I2(1) + (I2(2) - I2(1)) * rand()];
    disp(x_init);
%     x_init = [0.452077 0.648269243203752];


    y = x_init;
    v = [];

    for t = 1:T
        u = controller.getInputs(y(end,:));
        v = [v; u(1,:)];
        [~, x] = ode45(@engine_ode, [0 dt], y(end,:), ...
            odeset('abstol',1e-4,'reltol',1e-4), u(1,:));
        y = [y; x(end,:)];
    end

    % store trajectory
    trajectories{k} = y;
end

% eps = 0.000016; % approximate value of the ASF, obtained from the .ipynb file
% state_grid = SymbolicSet('engine_state_grid.bdd'); % SCOTS grid object
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
%             % Pick a random row (control input) from u_all
%             idx_rand = randi(size(u_all,1));
%             u_sel = u_all(idx_rand,:);
%         else
%             u_sel = controller.getInputs(current_y);
%             u_sel = u_sel(1,:);
%         end
% %         u_sel = controller.getInputs(current_y);
% %         u_sel = u_sel(1,:);
% 
%         v = [v; u_sel];
%         [~, x] = ode45(@engine_ode, [0 dt], current_y, ...
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

% plot target set
v = [0.4479 0.6473; 0.4559 0.6473; 0.4479 0.6553; 0.4559 0.6553];
patch('vertices', v, 'faces', [1 2 4 3], ...
      'facecolor', 'none', 'edgec', colors(2,:), 'linew', 1)
hold on

box on
grid on
axis([0.4479 0.4559 0.6473 0.6553])
axis tight
ylabel('$x_2$','Interpreter','latex', 'FontSize',18)
xlabel('$x_1$','Interpreter','latex', 'FontSize',18)

end

function dxdt = engine_ode(t, x, u)
    % Parameters
    a = 1.0/3.5;
    H = 0.18;
    l = 8;
    B = 2;
    W = 0.25; 
    psi = a + H*(1 + 1.5*(x(1)/W-1) - 0.5*(x(1)/W-1)^3);
       
    % Preallocate dxdt
    dxdt = zeros(size(x));

    % Dynamics
    dxdt(1) = 1/l*(psi - x(2)) + u(1);
    dxdt(2) = 1/(4*l*B*B)*(x(1) - u(2)*sqrt(x(2)));
end