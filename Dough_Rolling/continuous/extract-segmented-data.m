%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Segmented Data from tGau-BP-HMM Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
load('proc-data.mat')
% Remove first time-series (not used)
Xn_seg = {};
N = length(Xn);
for i = 2:N
    Xn_seg{i-1,1} = Xn{i};
end
N = length(Xn_seg);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Results from tGau-BP-HMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('seg-results.mat')
% Variables needed from tGau-BP-HMM Results
%   - Robotdata
%   - bestGauPsi
%   - bestGauPsiTrans
%   - groups

% Visualize Segmentation and Sigma-Clustering
close all; clc;
[ Segm_results Total_feats Clust_results Clust_feats my_color_map] = GetSegmentationResults(Robotdata, bestGauPsi, [1:Robotdata.N], 'Best estimated State Sequences', groups);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Segmented Trajectories with BP-HMM on 3D Cartesian Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select which sequences to visualize
seq = [1 3 5 7 10 12];

figure('Color',[1 1 1])

% Plot 3D Trajectories of Recordings with Colors Indicating Sequences
plotSegmentedData( Xn_seg, seq , Total_feats, Segm_results, my_color_map);

% Set Reference Frames
Robot_Base   = eye(4);
Rolling_Board = Table_Hn{1}(:,:,1);
visualizeRollingEnvironment(Robot_Base, Rolling_Board);

axis tight
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Recordings \n(Color Indicates Segments Extracted by BP-HMM)',length(seq)))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Segmented Trajectories with tGau-BP-HMM on 3D Cartesian Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select which sequences to visualize
seq = [1 3 5 7 10 12];

figure('Color',[1 1 1])

% Plot 3D Trajectories of Recordings with Colors Indicating Sequences
plotSegmentedData( Xn_seg, seq , Clust_feats, Clust_results, my_color_map);

% Set Reference Frames
Robot_Base   = eye(4);
Rolling_Board = Table_Hn{1}(:,:,1);
visualizeRollingEnvironment(Robot_Base, Rolling_Board);

axis tight
grid on 
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('%d Dough Rolling Recordings \n(Color Indicates Segments Extracted by tGau-BP-HMM)',length(seq)))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Most Likely Sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feats = length(Total_feats);
kappa = bestGauPsi.kappa;
A = zeros(feats,feats);
obs_seq = [];
for i=1:N
    tmp = bestGauPsi.Eta(i).eta;
    for ii=1:feats
        tmp(ii,ii) = tmp(ii,ii)/kappa;
    end
    A = A + tmp;
    obs_seq{i} = Segm_results{i}(:,1);
end

% Average Feature Transition Matrix
A = A/N;

% Compute Initial Distribution from Counts
% overall_patterns  = [];
possible_sequence = [];
counts = [];
k = 0;
for j=1:length(obs_seq)
    patterns = find_patterns(obs_seq{j});
    for i=1:length(patterns)
        if length(patterns{i,1}) == feats
            possible_sequence = [possible_sequence; patterns{i,1}];
            counts    = [counts; patterns{i,2}];
        end
    end
end
[unique_sequences, ids, seq_ids]=unique(possible_sequence,'rows');

for i=1:length(ids)
    seq_id_counts(i,1) = sum(counts(find(seq_ids==i)));
end
seq_ids_counts = seq_id_counts/sum(seq_id_counts);

% Conditional Probability of Possible Sequence Order
prob_seq = [];
for i=1:length(ids)
    poss_seq  = unique_sequences(i,:); 
    seq_prior = seq_ids_counts(i);
    seq_like = 1;
    for j=1:length(ids)-1
        seq_like  = seq_like + A(j,j+1) ; 
    end
    prob_seq(i,1) = seq_prior*seq_like;
end

[best_prob best_seq_id] = max(prob_seq);

% Task (Action) Sequence should be Executed in the following order:
task_sequence = unique_sequences(best_seq_id, :);

% Replace IDs with high-level actions (clusters)
tmp_task_sequence = [];
for i=1:length(task_sequence)
    feat =  task_sequence(i);    
    tmp_task_sequence = [tmp_task_sequence Clust_feats(feat)];
end

final_task_sequence = [];
for i=1:length(tmp_task_sequence)
    curr = tmp_task_sequence(i);
    if i == 1 || curr ~=tmp_task_sequence(i-1)
        final_task_sequence = [final_task_sequence curr]
    end    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Time-Series Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Per-Action Sequnces
action_sequences = [];
% for i=1:length(Xn_seg)
for i=1:1
    ts_seq = Clust_results{1};
    unique_seq = [];
    for i=2:length(ts_seq)
        if ts_seq(i-1,1)~=ts_seq(i,1)
            unique_seq = [unique_seq; ts_seq(i-1,:)];
        end
    end
    ind=strfind(reshape(unique_seq(:,1),1,[]),final_task_sequence);
    
    % Number of Repeated Sequences
    n_seq = length(ind);
    
    
end


% Per-Action


