%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Processing for Raw Data includes:
%       - Smoothing F/T Data
%       - Convert Orientation to Quaternions
%       - Fixing Discontinuities in Rotation Data
%       - Sub-Sampling data (Default: Raw Hz/5)
%       - Matching Hz for each different stream
%   Follow this script if you want to change the processed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
load('./mat/raw-data.mat')
data_dir = './mat/';
% Default processing parameters
sub_sampl_factor = 5;
N = length(data);

%% Data Pre-Processing
close all
N = length(data);
proc_data = {};

for jj=1:N
    
    % Pre-process data from passive arm
    options.match_hz    = 1;
    options.smooth_fact = 0.005;
    options.do_plot     = 1;
    options.tool_frame  = Hand_ft;
    options.base_frame  = [];
    options.title_name  = ['Sequence: ',data{jj}.name, ' ---- Arm: Passive'];
    [X_p, H_p, t_p] = preProcessEEData(data{jj}.passive.pose, data{jj}.passive.ft, ... 
        data{jj}.times, options);

    % Passive Arm Variables
    proc_data{jj}.passive.X = X_p;
    proc_data{jj}.passive.H = H_p;
    proc_data{jj}.passive.t = t_p;
    proc_data{jj}.passive.base_frame = eye(4);
    proc_data{jj}.passive.tool_frame = Hand_ft;
    
    % Pre-process data from active arm
    options.smooth_fact = 0.005;
    options.tool_frame  = Tool_ft;
    options.base_frame  = Active_robot;
    options.title_name  = ['Sequence: ',data{jj}.name, ' ---- Arm: Active'];
    [X_a, H_a, t_a] = preProcessEEData(data{jj}.active.pose, data{jj}.active.ft,  ...
        data{jj}.times, options);   
    
    % Active Arm Variables
    proc_data{jj}.active.X = X_a;
    proc_data{jj}.active.H = H_a;
    proc_data{jj}.active.t = t_a;    
    proc_data{jj}.active.base_frame = Active_robot;
    proc_data{jj}.active.tool_frame = Tool_ft;
    
    
    % Pre-process data from object features   
    options.smooth_fact = 0.01;   
    options.time_arm    = t_a;
    options.title_name  = ['Sequence: ',data{jj}.name, ' ----: Object'];    
    
    [RGB_feats] = preProcessObjectData(data{jj}.object.feats, options);   
    
    % Object Variables
    proc_data{jj}.object.feats = RGB_feats;
    
    % Global Task Variables
    proc_data{jj}.name         = data{jj}.name;
    proc_data{jj}.board_lframe = Cutting_board_l;
    proc_data{jj}.board_rframe = Cutting_board_r;
    
    % Plot Trajectories with Reference Frames etc
    plotBimanualTrajectories(proc_data{jj})
    
end

%% Save proc data to matfile
matfile = strcat(data_dir,'proc-data.mat');
save(matfile,'proc_data')

