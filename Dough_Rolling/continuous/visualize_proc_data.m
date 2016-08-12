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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Cartesian EE Trajectories in 3D Task Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Altogether with table and Reference Frames

figure('Color', [1 1 1])
for i=1:N
    cart_ee_traj = Xn{i};
    
    % Plot TrAjectories
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '-','Color',[rand rand rand],'LineWidth',1); hold on
    
    % Plot Starting and End Points
    start_points = [cart_ee_traj(1,1), cart_ee_traj(2,1),cart_ee_traj(3,1)];
    end_points   = [cart_ee_traj(1,end),cart_ee_traj(2,end),cart_ee_traj(3,end)];
    scatter3(start_points(:,1),start_points(:,2),start_points(:,3), 70, [0 1 0], 'filled'); hold on;    
    scatter3(end_points(:,1),end_points(:,2),end_points(:,3),70, [1 0 0], 'filled'); hold on;
    
    % Draw some frame of Start-end Trajectories
    H = EE_HnP{i};
    drawframe(H(:,:,1),0.02); 
    drawframe(H(:,:,end),0.02);
    
end


% Set Reference Frames
Robot_Base   = eye(4);
Rolling_Board = table_frame;
visualizeRollingEnvironment(Robot_Base, Rolling_Board);

axis tight
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Recordings (Color Indicates different time-series)',N))


% Per time-series
rc = ceil(sqrt(N));
figure('Color', [1 1 1])
for i=1:N
    subplot(rc,rc,i);
    cart_ee_traj = Xn{i};
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '-','Color',[rand rand rand],'LineWidth',1);
    grid on 
    xlabel('x');ylabel('y');zlabel('z');
    title(sprintf('Dough Rolling Recording %d',i))
end
suptitle('3-4 Rolling Sequences per Time-Series')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Full EE Data Time-Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Altogether
for i=1:N
    plotEEData(Xn{i}, [], sprintf('Dough Rolling Time Series %d',i)); 
end

