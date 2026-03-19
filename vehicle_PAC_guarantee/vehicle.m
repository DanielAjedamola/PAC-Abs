%
% vehicle.m
%
% created on: 24.02.2025
%     author: Daniel A.
% see readme file for more information on the vehicle example
%
% you need to run ./vehicle binary first 
%
% so that the files: vehicle_ss.bdd 
%                    vehicle_obst.bdd
%                    vehicle_target.bdd
%                    vehicle_controller.bdd 
% are created
%

function vehicle
clear set
close all



%% simulation

N = 100; % Number of trajectories
tau = 0.25;

% Load controller and target region
controller = SymbolicSet('vehicle_controller.bdd', 'projection', [1 2 3]);
target = SymbolicSet('vehicle_target.bdd');

% Initialize storage for all trajectories
Y = cell(N, 1); % Stores trajectories
V = cell(N, 1); % Stores control inputs

for n = 1:N
    % Sample initial state x0 randomly from [0, 0.2]^3
    x0 = [rand * 0.2, rand * 0.2, rand * 0.2];
%     x0 = [0, 0.2, 0];
    
    % Initialize trajectory and control input storage
    y = x0;
    v = [];

    % Generate trajectory using Forward Euler
    while (1)
        if target.isElement(y(end, :))
            break;
        end 

        u = controller.getInputs(y(end, :));
        v = [v; u(1, :)];

        % Forward Euler Integration
        xk = y(end, :).';  % Convert row vector to column
        dxdt = unicycle_ode(0, xk, u(1, :));
        x_next = xk + tau * dxdt;  % Euler step

        y = [y; x_next.'];  % Store the new state as a row vector
    end

    % Store trajectory and control inputs
    Y{n} = y;
    V{n} = v;
end

disp('All trajectories complete');

%% Plot the vehicle domain
colors = get(groot, 'DefaultAxesColorOrder');
figure('DefaultAxesFontSize', 16)
figure(1)

% Load and plot the symbolic set containing the abstract state space
set = SymbolicSet('vehicle_ss.bdd', 'projection', [1 2]);
plotCells(set, 'facecolor', 'none', 'edgec', [0.8 0.8 0.8], 'linew', .1)
hold on

% Load and plot the symbolic set containing obstacles
set = SymbolicSet('vehicle_obst.bdd', 'projection', [1 2]);
plotCells(set, 'facecolor', colors(1, :) * 0.5 + 0.5, 'edgec', colors(1, :), 'linew', .1)

% Plot the real obstacles and target set
plot_domain

% Load and plot the symbolic set containing the target set
set = SymbolicSet('vehicle_target.bdd', 'projection', [1 2]);
plotCells(set, 'facecolor', colors(2, :) * 0.5 + 0.5, 'edgec', colors(2, :), 'linew', .1)

% Generate a set of N distinct colors
traj_colors = lines(N);

% Plot all trajectories
for n = 1:N
    y = Y{n};
    plot(y(:, 1), y(:, 2), '.-', 'Color', traj_colors(n, :), 'LineWidth', 1.5) % Plot trajectory
    plot(y(1, 1), y(1, 2), '.', 'Color', traj_colors(n, :), 'color', colors(5, :), 'markersize', 20) % Plot initial state
end

box on
axis([-.5 10.5 -.5 10.5])
axis tight
ylabel('$y(m)$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('$x(m)$', 'Interpreter', 'latex', 'FontSize', 18)
end


function dxdt = unicycle_ode(t,x,u)
  dxdt = zeros(3,1);
  c=atan(tan(u(2))/2);
  dxdt(1)=u(1)*cos(c+x(3))/cos(c);
  dxdt(2)=u(1)*sin(c+x(3))/cos(c);
  dxdt(3)=u(1)*tan(u(2));
end


function plot_domain

colors=get(groot,'DefaultAxesColorOrder');

v=[9 0; 9.5  0; 9 0.5; 9.5 .5];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(2,:),'edgec',colors(2,:));


v=[1     0  ;1.2  0   ; 1     9    ; 1.2 9   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2   0  ;2.4  0   ; 2.2   5    ; 2.4 5   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2   6  ;2.4  6   ; 2.2   10   ; 2.4 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[3.4   0  ;3.6  0   ; 3.4   9    ; 3.6 9   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[4.6   1  ;4.8  1   ; 4.6   10   ; 4.8 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.8   0  ;6    0   ; 5.8   6    ; 6   6   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.8   7  ;6    7   ; 5.8   10   ; 6   10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[7     1  ;7.2  1   ; 7     10   ; 7.2 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.2   0  ;8.4  0   ; 8.2   8.5  ; 8.4 8.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   8.3;9.3  8.3 ; 8.4   8.5  ; 9.3 8.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   7.1;10   7.1 ; 9.3   7.3  ; 10  7.3 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   5.9;9.3  5.9 ; 8.4   6.1  ; 9.3 6.1 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   4.7;10   4.7 ; 9.3   4.9  ; 10  4.9 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   3.5;9.3  3.5 ; 8.4   3.7  ; 9.3 3.7 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   2.3;10   2.3 ; 9.3   2.5  ; 10  2.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


end