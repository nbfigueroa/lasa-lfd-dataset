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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Cartesian EE Trajectories in 3D Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% Altogether
figure('Color', [1 1 1])
for i=1:N
    
    cart_ee_traj = EE_CART_Pn{i};
    % Plot TrAjectories
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '-','Color',[rand rand rand],'LineWidth',1); hold on
     
    % Plot Starting and End Points
    start_points = [cart_ee_traj(1,1), cart_ee_traj(2,1),cart_ee_traj(3,1)];
    end_points   = [cart_ee_traj(1,end),cart_ee_traj(2,end),cart_ee_traj(3,end)];
    scatter3(start_points(:,1),start_points(:,2),start_points(:,3), 70, [0 1 0], 'filled'); hold on;    
    scatter3(end_points(:,1),end_points(:,2),end_points(:,3),70, [1 0 0], 'filled'); hold on;
    
    % Draw some frame of Start-end Trajectories
    H = EE_Hn{i};
    drawframe(H(:,:,1),0.02); 
    drawframe(H(:,:,end),0.02);
     
end
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Recordings (Color Indicates different time-series)',N))

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
    cart_ee_traj = EE_CART_Pn{i};
    plot3(cart_ee_traj(1,:), cart_ee_traj(2,:),cart_ee_traj(3,:), '-','Color',[rand rand rand],'LineWidth',1);
    grid on 
    xlabel('x');ylabel('y');zlabel('z');
    title(sprintf('Dough Rolling Recording %d',i))
end
suptitle('3-4 Rolling Sequences per Time-Series')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Full EE Data Time-Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
% All
for i=1:N
    X = [EE_CART_Pn{i};EE_CART_On{i};EE_FTn{i}];
    plotEEData(X, [], sprintf('Dough Rolling Time Series %d',i)); 
end
