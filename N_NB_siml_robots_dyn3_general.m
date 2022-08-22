% This script is for N noelbots (DDWMR vehicle paltoon)...
...simulink model and hardware robot run file

% Run_this_file_for_N_nolebots_vehicle_following

clear all
close all
clc

% % % %% master ip address for connection
% rosshutdown
% % rosinit('192.168.2.14')
% rosinit

%% parameter for robot model
r  = 0.05;                       % wheel's radium (m)
mc = 10;                         % robot platform's weight (kg)
mw = 2.5;                        % Wheel's weight (kg)
m  = mc+2*mw;                    % total weight (kg)
L  = 0.235/2;                    % half of distance between two wheels (m)
d  = 0.1;                        % distance betwen medium point of axis of wheels and center of mass(m)
Ic = 0.05;                       % inertia of robot about vertical axis (kg.m^2)
Iw = 0.025;                      % inertia of wheel about wheel axis (kg.m^2)
Im = 0.025;                      % inertia of wheel about wheel diameter (kg.m^2)
I  = Ic +mc*d^2+2*mw*L^2+2*Im;   % inertia of whole robot (kg.m^2)

%% Vehicle following paramerters

% Controller parameters
% h = 0.5;
h = 0.7;   % robot headway
tau = 0.1; % time constant 

% mobile robot gain values
kp = 0.2;   % proportional gain
kd1 = 0.7;  % differrntial gain
kd2 = 0;
vehicle_length_1 = 0.35;    % robot length
stand_still_1 = 0.65;       % stand still distance 

% for vehicle following
dim = 2; % Number of states per agent (This simulation contains homogenous agents)
n  = 10; % total number of nolebots
nr = 3; % number of real robots (Nolebots)
ns = n-nr; % number of following nolebot simulink models 
r_ind = zeros(1,n); % initializing the array to store boolean values based on the robots position index
LR = (vehicle_length_1 + stand_still_1)*ones(n-1,1); % length of robot model and desired distance b/w robots
% r_LR = (vehicle_length_1 + stand_still_1)*ones(nr-1,1); % % length of robot model and disered distance b/w robots (not considering first robot)

%% indicator based on robots position

fprintf("Give position index of %d robots in %d vehicle string \n", nr, n);
% storing the user input position indices of the nr number of robots
robot_pos_index = [1,6,10];
% creating a boolean array, considering 0 at robot position and 1 at
% virtual robot position
for i=1:n
    % putting 0 or 1 at robot position and 1 at virtual robot position
    if any(robot_pos_index(:) == i)
        r_ind(i) = 0;
    else
        r_ind(i) = 1;
    end
end

isko = 1:numel(r_ind);  % to store the indeces of total number of nolebots (both model and robot)

r_ind_for_A1 = r_ind(2:end);
r_ind_for_A2 = r_ind(1:end-1);
r_ind_for_A3 = r_ind(2:end-1);


r_ind_b1 = zeros(n-1,numel(robot_pos_index)-1);
r_ind_b2 = zeros(n-1,numel(robot_pos_index));

if robot_pos_index(1) == 1
    for i = 2:numel(robot_pos_index)
        r_ind_b1((robot_pos_index(i)-1),i-1) = 1;
    end
else
     for i = 1:numel(robot_pos_index)
        r_ind_b2((robot_pos_index(i)-1),i) = 1;
     end
end

if robot_pos_index(1) == 1
    r_ind_b = [zeros(n-1,1), r_ind_b1];    
else
    r_ind_b = r_ind_b2;    
end

ur_choose_idx = find(1-r_ind);

%% Initial condition

% initial conditions of the robot model
theta_0=0;
x_0=0;  % initial x position of model
y_0=0;  % initial y position of model
v_0=0;  % initial x velocity of model
w_0=0;  % initial angular velocity about z of model

% initialization of states (distance between robots and velocities of robots)
d_init  = 2*LR(1)*linspace(1,n-2,n-1).'; % initial condition of distance between all 
v_init  = zeros(n-1,1);   % initial condition of velocity of robots
x_init  = kron(d_init,[1;0]) + kron(v_init,[0;1]); % intital states of simulink model
xc_init = zeros(n-1,1);  % initial conditions of controller states

if robot_pos_index(1) == 1
    s_R1_init = 0;
else
    s_R1_init = -sum(d_init(1:(robot_pos_index(1)-1)));
end
s_R2_init = -sum(d_init(1:(robot_pos_index(2)-1)));  % robot-1 initial position 
s_R3_init = -sum(d_init(1:(robot_pos_index(3)-1)));  % robot-2 initial position 

% inital conditions of 3 robots (need to change for general case)
q_initial_1 = [theta_0; s_R1_init; y_0];   % initial for reference kinematic model
V_initial_1 = [v_init(1);w_0];
q_initial_2 = [theta_0; s_R2_init; y_0];   % initial for reference kinematic model
V_initial_2 = [v_init(2);w_0];
q_initial_3 = [theta_0; s_R3_init; y_0];   % initial for reference kinematic model
V_initial_3 = [v_init(3);w_0];


%% Generating matrices

% creating a matrix A 
A = (1/tau)*(kron(diag(r_ind_for_A1)*diag(ones(n-1,1)),[0 -tau; 0 -1])...
    + kron(diag((r_ind_for_A3),-1),[0 tau;0 0]));
% coefficient matrix of the velocity vector (real robots or models) 
B_v_init = kron(diag((1 - r_ind_for_A1)),[-tau; -1]) + kron(diag((1 - r_ind_for_A3),-1),[tau; 0]);
idx_nonzerolines = sum(abs(B_v_init),1)>0;
B_v_leading = zeros(2*(n-1),1);
B_v_leading(1) = tau;
B_v = (1/tau)*[B_v_leading, B_v_init(:,idx_nonzerolines)];
% coefficient matrix of the control input vector (model)
B_u = (1/tau)*kron(eye(n-1),[0; 1]);
% C matrix for output
C = eye((n-1)*dim);   
D = zeros(size(C,1),size(B_u,2));
% Transforms from [d1 v1 d2 v2 ... dn vn]' to [d1 d2 ... dn v1 v2 ... vn];
C_transform = [kron(eye(n-1),[1 0]);  
               kron(eye(n-1),[0 1])]; 
Ce          = [eye(n-1) -h*eye(n-1)]*C_transform; 

%% controller

% A matrix for controller dynamics (any number of models and robot at any position)
A_con = (1/h)*(-eye(n-1) + diag(ones(n-2,1),-1));
% A_con = (1/h)*(-eye(n-1) + diag(ones(n-2,1),-1) - 0.2*diag((1-r_ind(2:end)))*diag(ones(n-2,1),-1));
B_con = zeros((n-1),1);
B_con(1)=(1/h);

C_transform_rev = [kron(eye(n-1),[1; 0])    kron(eye(n-1),[0; 1])]; 
