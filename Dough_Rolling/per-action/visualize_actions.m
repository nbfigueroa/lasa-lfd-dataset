%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Dough-Rolling actions contain end-effector:
% - positions:         actions.data{i,j}(1:3,:)   (3-d: x, y, z)
% - orientations:      actions.data{i,j}(4:7,:)   (4-d: q_i, q_j, q_k, q_w)
% - forces:            actions.data{i,j}(8:10,:)  (3-d: f_x, f_y, f_z)
% - torques:           actions.data{i,j}(11:13,:) (3-d: tau_x, tau_y, tau_z)
%
%   i: action instance, j:
%
% PROC DATA Sampling Rate: 100Hz
%   
% Details:
%    - Each time-series contains 3 - 4 iterations of a rolling sequence, 
%      consisting of a sequence of sub-actions:
%      1.- reach
%      2.- roll
%      3.- back
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;
load('actions.mat')
action_sequence = actions.sequence;
action_data = actions.data;
N = length(action_data);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize each action in 3D Cartesian Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color', [1 1 1])
for j=1:length(action_sequence)
    subplot(1,length(action_sequence), j)
    action = j;
    
    for i=1:N
        action_ts = action_data{i,action};
        plot3(action_ts(1,:),action_ts(2,:),action_ts(3,:),'-o','MarkerSize',0.5);hold on;
    end

    % Set Reference Frames
    Robot_Base   = eye(4);
    Rolling_Board = table_frame;
    visualizeRollingEnvironment(Robot_Base, Rolling_Board);
    grid on
    axis equal
    xlabel('x');ylabel('y');zlabel('z');
    title(sprintf('%d time-series for Action id: %d',N,action_sequence(action)))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize each action as EE Time Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotEEActionInstances(action_data, action_sequence)
