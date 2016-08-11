%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Processing for Raw Data includes:
%       - Smoothing F/T Data
%       - Convert Orientation to Quaternions
%       - Fixing Discontinuities in Rotation Data
%       - Sub-Sampling data (Default: Raw Hz/5)
%   Follow this script if you want to change the processed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
load('raw-data.mat')

% Default processing parameters
sub_sampl_factor = 5;
N = length(EE_CART_Pn);

% Type of Rotation
[ori, ~] = size(EE_CART_Pn{1});

%% %%%%%%%%%%%%%%%%%%%%%%
% Sub Sample RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%

EE_PnP = {};
EE_OnP = {};
EE_FTnP = {};
EE_HnP = {};

for jj=1:N
    P_tmp = EE_CART_Pn{jj,1};
    P_sample = P_tmp(:,1:sub_sampl_factor:end);
    
    O_tmp = EE_CART_On{jj,1};
    O_sample = O_tmp(:,1:sub_sampl_factor:end);
    
    FT_tmp = EE_FTn{jj,1};
    FT_sample = FT_tmp(:,1:sub_sampl_factor:end);
    
    H_tmp = EE_Hn{jj,1};
    H_sample = H_tmp(:,:,1:sub_sampl_factor:end);
    
    EE_PnP{jj,1} = P_sample;
    EE_OnP{jj,1} = O_sample;
    EE_FTnP{jj,1} = FT_sample;
    EE_HnP{jj,1} = H_sample; 
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concat data to 13-d time-series (X) and process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xn = [];
for ii=1:N
    
    % Extract quaternion and check discont
    H = EE_HnP{ii,1};
    R = H(1:3,1:3,:);
    % Check Rotations twice!
    Q = checkRotations(quaternion(R));
    EE_Q = checkRotations(Q);
    
    % Concat Pos + Orient
    X = cat(1,EE_PnP{ii,1},EE_Q);    
    
    %From ee_ft
    FT = cat(1,EE_FTnP{ii,1}(1:3,:),EE_FTnP{ii,1}(4:6,:));      
    
    % Smoothen FT Sensor Data
    sm_ft = FT;
    for kk=1:size(sm_ft,1)
        sm_ft(kk,:) = smooth(sm_ft(kk,:),0.01,'moving');
    end
    X_tmp = cat(1,X,sm_ft);  
    Xn{ii,1} = X_tmp;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize Full EE Data Time-Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
for i=1:N
    plotEEData(Xn{i}, [], sprintf('Dough Rolling Time Series %d',i)); 
end
