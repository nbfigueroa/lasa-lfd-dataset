%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Dough-Rolling sequence contain end-effector:
% - positions:         action_sequence{i}(1:3,:)   (3-d: x, y, z)
% - orientations:      action_sequence{i}(4:7,:)   (4-d: q_i, q_j, q_k, q_w)
% - forces:            action_sequence{i}(8:10,:)  (3-d: f_x, f_y, f_z)
% - torques:           action_sequence{i}(11:13,:) (3-d: tau_x, tau_y, tau_z)
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
load('action-sequences.mat')
N = length(action_sequences);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Extracted Rolling Sequences from tGau-BP-HMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the Sequences to Visualize
M = length(action_sequences);
seq = [1:5:M];

%%%%%%%%%%  Visualize 3D trajectories %%%%%%%%%%
figure('Color',[1 1 1])
clear actions;
for i=1:length(seq)
  sequence = action_sequences{seq(i)};  
  for j=1:length(inf_action_sequence)
    action            = inf_action_sequence(j);
    action_ids        = find(sequence(size(sequence,1),:) == action);
    action_data       = sequence(1:end-1,action_ids);    
    plot3(action_data(1,:),action_data(2,:),action_data(3,:), 'Color', my_color_map(action,:)); hold on;
  end
end
% Set Reference Frames
Robot_Base   = eye(4);
Rolling_Board = table_frame;
visualizeRollingEnvironment(Robot_Base, Rolling_Board);

axis tight
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Sequences extracted with tGau-BP-HMM \n (Color Indicates different actions)',length(seq)))


%% %%%%%%%%  Visualize Full EE Variables as Time_Series %%%%%%%%%%
seq = [1:5:M];
for i=1:length(seq)
    plotLabeledEEData(action_sequences{seq(i)}, [], sprintf('Dough Rolling Sequence %d',seq(i))); 
end
