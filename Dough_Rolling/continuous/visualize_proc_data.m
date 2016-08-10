%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Dough-Rolling processed (smoothed f/t trajactories, fixed rotation
% discontinuities and sub-sampled) data contains end-effector:
% - positions:         Xn{i}(1:3,:)   (3-d: x, y, z)
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
load('proc-data.mat')
N = length(Xn);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Cartesian EE Trajectories in 3D Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Altogether
figure('Color', [1 1 1])
for i=1:N
    cart_ee_traj = Xn{i};
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '--','Color',[rand rand rand],'LineWidth',2); hold on
end
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Recordings (Color Indicates different time-series)',N))

% Per time-series
rc = ceil(sqrt(N));
figure('Color', [1 1 1])
for i=1:N
    subplot(rc,rc,i);
    cart_ee_traj = Xn{i};
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '--','Color',[rand rand rand],'LineWidth',2);
    grid on 
    xlabel('x');ylabel('y');zlabel('z');
    title(sprintf('Dough Rolling Recording %d',i))
end
suptitle('3-4 Rolling Sequences per Time-Series')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Full EE Data Time-Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Altogether
for i=1:N
    plotEEData(Xn{i}, [], sprintf('Dough Rolling Time Series %d',i)); 
end
