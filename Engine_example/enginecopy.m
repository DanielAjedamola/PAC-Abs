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
% clearvars -except hFig  % keeps figure handle if it exists
% set(0, 'DefaultFigureWindowStyle', 'normal'); % ensure figure is visible
% 
% % Create or reuse existing figure
% if ~exist('hFig', 'var') || ~isvalid(hFig)
%     hFig = figure('Name', 'Engine Trajectories', 'NumberTitle', 'off');
% else
%     figure(hFig);
%     hold on;  % keep adding to existing plot
% end

%% Simulation setup
I1 = [0.4479, 0.4540];   % interval for x(1)
I2 = [0.6474, 0.6553];   % interval for x(2)
N = 1;                   % number of trajectories
controller = SymbolicSet('engine_controller.bdd');
T = 30000; 
dt = 0.1;

%% Plot setup
colors = get(groot,'DefaultAxesColorOrder');

% Plot controller domain (gray background)
x = controller.points;
plot(x(:,1), x(:,2), '.', 'color', [0.8 0.8 0.8])
hold on
box on; grid on;
axis([0.4475 0.4565 0.6469 0.6556])
ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 18)

% Plot target set
v = [0.4479 0.6473; 0.4559 0.6473; 0.4479 0.6553; 0.4559 0.6553];
patch('vertices', v, 'faces', [1 2 4 3], ...
      'facecolor', 'none', 'edgec', colors(2,:), 'linew', 1)

%% Sequential trajectory simulation and plotting
for k = 1:N
    % Random initial condition from I1 × I2
    x_init = [I1(1) + (I1(2) - I1(1)) * rand(), ...
              I2(1) + (I2(2) - I2(1)) * rand()];
%     disp(['Initial condition ', num2str(k), ': ', num2str(x_init)])

    y = x_init;
    v = [];

    % Simulate one trajectory
    for t = 1:T
        u = controller.getInputs(y(end,:));
        v = [v; u(1,:)];
        [~, x_next] = ode45(@engine_ode, [0 dt], y(end,:), ...
            odeset('abstol', 1e-4, 'reltol', 1e-4), u(1,:));
        y = [y; x_next(end,:)];
    end

    % Plot trajectory immediately
    plot(y(:,1), y(:,2), '.-', 'color', colors(mod(k-1, size(colors,1))+1,:), 'markersize', 0.5)
    % plot(y(1,1), y(1,2), '.', 'color', colors(5,:), 'markersize', 15)
    plot(y(1,1), y(1,2), 'k.', 'markersize', 15)
    drawnow  % update plot incrementally
end

% %% Sequential trajectory simulation and plotting (with perturbation and ASF)
% eps = 0.000016;  % approximate ASF value
% state_grid = SymbolicSet('engine_state_grid.bdd');  % SCOTS grid object
% abs_states = state_grid.points;                     % all abstract states
% 
% for k = 1:N
%     % Initial condition: perturbation around x0
%     x_init = [I1(1) + (I1(2) - I1(1)) * rand(), ...
%               I2(1) + (I2(2) - I2(1)) * rand()];
%     y = x_init;
%     v = [];
% 
%     % Simulate one trajectory
%     for t = 1:T
%         current_y = y(end,:);
% 
%         % Find nearby abstract states
%         idx = vecnorm(abs_states - current_y, 2, 2) <= eps;
%         u_all = [];
%         for i = find(idx)'
%             u_tmp = controller.getInputs(abs_states(i,:));
%             u_all = [u_all; u_tmp(1,:)];
%         end
% 
%         % Select control input
%         if ~isempty(u_all)
%             u_sel = u_all(randi(size(u_all,1)), :);  % random among possible
%         else
%             u_sel = controller.getInputs(current_y);
%             u_sel = u_sel(1,:);
%         end
% 
%         v = [v; u_sel];
% 
%         % Simulate dynamics
%         [~, x_next] = ode45(@engine_ode, [0 dt], current_y, ...
%             odeset('abstol', 1e-4, 'reltol', 1e-4), u_sel);
%         y = [y; x_next(end,:)];
%     end
% 
%     % Store trajectory (optional)
%     trajectories{k} = y;
% 
%     % Plot trajectory
%     plot(y(:,1), y(:,2), '.-', 'color', colors(mod(k-1, size(colors,1))+1,:), 'markersize', 0.5)
%     hold on
%     plot(y(1,1), y(1,2), 'k.', 'markersize', 15)   % black initial point
%     drawnow
% end

end

%% Engine dynamics
function dxdt = engine_ode(t, x, u)
    a = 1.0/3.5;
    H = 0.18;
    l = 8;
    B = 2;
    W = 0.25; 
    psi = a + H*(1 + 1.5*(x(1)/W - 1) - 0.5*(x(1)/W - 1)^3);

    dxdt = zeros(size(x));
    dxdt(1) = 1/l*(psi - x(2)) + u(1);
    dxdt(2) = 1/(4*l*B*B)*(x(1) - u(2)*sqrt(x(2)));
end