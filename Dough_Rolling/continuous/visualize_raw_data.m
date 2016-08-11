%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Dough-Rolling raw data contains end-effector:
% - positions:         EE_CART_Pn (3-d: x, y, z)
% - orientations:      EE_CART_On (3-d: roll, pitch, yaw)
% - forces/torques:    EE_FTn     (6-d: f_x, f_y, f_z, tau_x, tau_y, tau_z)
% - full 6-dof matrix: EE_Hn      (4x4: R \in R^{3x3}, t )
% We also include the relative frame of the table wrt. the robot base. (Table_Hn)
%
% RAW DATA Sampling Rate: 500Hz
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
load('raw-data.mat')
N = length(EE_CART_Pn);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Cartesian EE Trajectories in 3D Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Altogether
figure('Color', [1 1 1])
for i=1:N
    cart_ee_traj = EE_CART_Pn{i};
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '-*','Color',[rand rand rand],'LineWidth',0.5,'MarkerSize',1); hold on
end
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Recordings (Color Indicates different time-series)',N))

% Per time-series
rc = ceil(sqrt(N));
figure('Color', [1 1 1])
for i=1:N
    subplot(rc,rc,i);
    cart_ee_traj = EE_CART_Pn{i};
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '-*','Color',[rand rand rand],'LineWidth',0.5,'MarkerSize',1);
    grid on 
    xlabel('x');ylabel('y');zlabel('z');
    title(sprintf('Dough Rolling Recording %d',i))
end
suptitle('3-4 Rolling Sequences per Time-Series')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Full EE Data Time-Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
% Altogether
for i=1:N
    X = [EE_CART_Pn{i};EE_CART_On{i};EE_FTn{i}];
    plotEEData(X, [], sprintf('Dough Rolling Time Series %d',i)); 
end
