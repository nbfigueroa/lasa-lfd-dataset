%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order to learn a probabilistic representation of our actions
% we must first process the data a bit. 
% We use the data directly extracted from the tGau-BP-HMM.
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
clear all;clc;close all;
load('actions.mat')
action_sequence = actions.sequence;
action_data = actions.data;
N = length(action_data);
%% %%%%%%%%%%%%%%%%%%%%%
% Pre-processing station
%%%%%%%%%%%%%%%%%%%%%%%%
% TODO...
action = 1;

% Check for time-series time lengths
action_lengths = zeros(N,1);
for i=1:N
    action_lengths(i) = length(action_data{i,action});
end

figure('Color',[1 1 1])
hist(action_lengths);hold on
mode(action_lengths)
std(action_lengths)

% Disclaimer: The following pre-processing steps are considering that we 
% will learn a Dynamical System for pos/orientation encoded as a GMM.
% The most important preprocessing step is to check that all trajectories 
% end at 0 velocity; i.e. converge at the attractor.
min_lengths = length(Phase1(1).in_dough.EE_POS);
for ii=1:length(Phase1)   
        min_lengths = [min_lengths length(Phase1(ii).in_dough.EE_POS)];
end 
mean_length = ceil(mean(min_lengths));
min_length = min(min_lengths);
k = 1;
for jj=1:9
    figure('Color',[1 1 1])
    for ii=1:5
        subplot(5,1,ii)
        Phase = Phase1(k); 
        POS = Phase.in_dough.EE_POS;
        New_POS = interpolateSpline(POS,min_length);

        ORI = Phase.in_dough.EE_ORI;
        New_ORI = interpolateSpline(ORI,min_length);
        
        FOR = Phase.in_dough.EE_FOR;
        New_FOR = interpolateSpline(FOR,min_length);
        
        TQS = Phase.in_dough.EE_TQS;
        New_TQS = interpolateSpline(TQS,min_length);
        
        Interp_Data = [New_POS; New_ORI; New_FOR; New_TQS];
        plot(Interp_Data')
        name=strcat('Sequence',num2str(k))
        title(name)
        Phase1(k).in_dough_interp.EE_POS = New_POS;
        Phase1(k).in_dough_interp.EE_ORI = New_ORI;
        Phase1(k).in_dough_interp.EE_FOR = New_FOR;
        Phase1(k).in_dough_interp.EE_TQS = New_TQS;
        k=k+1;
    end
end

