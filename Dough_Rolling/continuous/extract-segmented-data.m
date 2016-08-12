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
possible_sequence = [];
counts = [];
k = 0;
for j=1:length(obs_seq)
    patterns = find_patterns(obs_seq{j});
    for i=1:length(patterns)
        if length(patterns{i,1}) == feats % Could be >2 to find sub-sequences
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

% Task (Feat) Sequence should be Executed in the following order:
task_sequence = unique_sequences(best_seq_id, :);

% Replace Feat IDs with high-level Actions (clusters)
tmp_task_sequence = [];
for i=1:length(task_sequence)
    feat =  task_sequence(i);    
    tmp_task_sequence = [tmp_task_sequence Clust_feats(feat)];
end

final_task_sequence = [];
for i=1:length(tmp_task_sequence)
    curr = tmp_task_sequence(i);
    if i == 1 || curr ~=tmp_task_sequence(i-1)
        final_task_sequence = [final_task_sequence curr];
    end    
end

display(final_task_sequence)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Time-Series Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; close all;

% Unique Cluster Assignments per Time-Series
uni_clust_results = [];
for ii=1:length(Xn_seg)
    ts_seq = Clust_results{ii};
    unique_seq = [];
    for i=2:length(ts_seq)
        if ts_seq(i-1,1)~=ts_seq(i,1)
            unique_seq = [unique_seq; ts_seq(i-1,:)];
        end
    end    
    uni_clust_results{ii} = unique_seq;    
end


% Per-task sequences (Same action sequence in each time-series)
inf_action_sequence = final_task_sequence;
Xn_seg = Xn_seg;

% extractActionSequences(inf_action_sequence, uni_cluster_results, Xn_seg)
action_sequences = [];
k = 0;
for ii=1:length(Xn_seg)
    unsegmented_ts = Xn_seg{ii};
    ts_action_assign = Uni_clust_results{ii};
    
    % Find the ordered sequence withing time-series
    ind=strfind(reshape(ts_action_assign(:,1),1,[]),inf_action_sequence);
    
    % Number of Repeated Sequences
    n_seq = length(ind);
    n_act = length(inf_action_sequence);
        
    % Extract Sequences    
    for jj=1:n_seq
        k = k + 1;
        if ind(jj) == 1        
            start_seq = 1            
        else
            start_seq =  ts_action_assign(ind(jj)-1,2)
        end
        end_seq   = ts_action_assign(ind(jj)+n_act-1,2)        
                
        act_assign = ts_action_assign(ind(jj):ind(jj)+n_act-1,:);
        one_sequence_in_ts = unsegmented_ts(:,start_seq:end_seq);
        seq_offset = act_assign(end,2) - length(one_sequence_in_ts)
        act_assign(:,2) = act_assign(:,2) - seq_offset;
        act_segms = [];
        for kk=1:length(inf_action_sequence)
            if kk==1
                start_action = 1;
            else
                start_action = act_assign(kk-1,2)
            end
            act_segms = [act_segms; start_action act_assign(kk,2)]
        end
        
        act_labels = [];
        act_labels = [act_labels inf_action_sequence(1)];
        for kk=1:length(inf_action_sequence)
            act_labels = [act_labels ones(1,act_segms(kk,2)-act_segms(kk,1))*inf_action_sequence(kk)];
        end
        one_sequence_in_ts = [one_sequence_in_ts;act_labels];
        action_sequences{k,1} = one_sequence_in_ts;        
    end
end


figure('Color',[1 1 1])
rc = ceil(sqrt(length(Xn_seg)));
plot3(one_sequence_in_ts(1,:),one_sequence_in_ts(2,:),one_sequence_in_ts(3,:),'-.','LineWidth',2); hold on
xlabel('x');ylabel('y');zlabel('z');
title(sprintf('Time-Series %d',ii))
grid on;
axis tight;

% Per-Action (each time series is an individual action)
% extractActions(inf_action_sequence, uni_cluster_results, Xn_seg)

