%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Dough-Rolling sequence contain end-effector:
% - positions:         act{i}(1:3,:)   (3-d: x, y, z)
% - orientations:      Xn{i}(4:6,:)   (3-d: roll, pitch, yaw)
% - forces:            Xn{i}(7:9,:)   (3-d: f_x, f_y, f_z)
% - torques:           Xn{i}(10:12,:) (3-d: tau_x, tau_y, tau_z)
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

clc;
% Select the Sequences to Visualize
M = length(action_sequences);
seq = [1:5];

%%%%%%%%%%  Visualize 3D trajectories %%%%%%%%%%
figure('Color',[1 1 1])
for i=1:length(seq)
  sequence = action_sequences{seq(i)};  
  plot3(sequence(1,:),sequence(2,:),sequence(3,:)); hold on;
end

% Set Reference Frames
Robot_Base   = eye(4);
Rolling_Board = Table_Hn{1}(:,:,1);
visualizeRollingEnvironment(Robot_Base, Rolling_Board);

axis tight
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Sequences extracted with tGau-BP-HMM \n (Color Indicates different sequences)',length(seq)))


%%%%%%%%%%  Visualize Full EE Variables as Time_Series %%%%%%%%%%

for i=1:length(seq)
    plotEEData(action_sequences{seq(i)}, [], sprintf('Dough Rolling Sequence %d',seq(i))); 
end
