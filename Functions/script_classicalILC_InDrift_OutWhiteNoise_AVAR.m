%%%%%%%%%%%%% script_classicalILC_InDrift_OutWhiteNoise_AVAR.m %%%%%%%%%%%%
%% Purpose:
%   To demonstrate ILC using the example from 
%   K. L. Moore, Y. Chen and H. -S. Ahn, "Iterative Learning Control: 
%   A Tutorial and Big Picture View," Proceedings of the 45th IEEE 
%   Conference on Decision and Control, San Diego, CA, USA, 2006.
%
%   Reference Input: Sinusoid/PRBS + RW in iteration domain
%   Measurement Noise: White Noise in sample domain
%
% Author:  Satya Prasad
% Updated: 2025/02/10

%% Prepare the workspace
clearvars
close all
clc

%% Define inputs and other parameters
number_of_iterations = 2^16;
number_of_runs = 1;
flag_input = 2; % 1 -> Sinusoid input, 2 -> PRBS input

sampling_frequency  = 50; % Sampling frequency [Hz]
sampling_interval   = 1/sampling_frequency; % Sampling interval/Simulation step size [s]
start_time          = 0; % Seconds
end_time            = 100*(2*pi*(100/8)*sampling_interval); % Seconds
time_vector         = start_time:sampling_interval:end_time;
length_of_iteration = length(time_vector); % Length for which one iteration lasts
dis_time_vector     = 0:length_of_iteration-1;

%% Noise parameters
power_spectral_density  = 4e-4; % [unit^2 s]
random_walk_coefficient = 0.02; % [unit/sqrt(s)] 0.2, 0.1, 0.02, 0.01
file_name = "20250210_InDrift0p02_OutWhite0p0004_AVAR.mat";

%% Plant definition
A = [-0.8, -0.22;...
        1, 0];
B = [0.5, 1]';
C = [1, 0.5];
D = 0;

%% Desired trajectory
if 1==flag_input
    reference_input = sin(8/100*dis_time_vector)'; % Desired trajectory to be tracked
elseif 2==flag_input
    rng('default')
    coin_toss   = rand(length_of_iteration,1)';
    flip_sign   = coin_toss>.95;
    total_flips = cumsum(flip_sign);
    reference_input = mod(total_flips,2)'; % Desired trajectory to be tracked
end % NOTE: END IF statement 'flag_input'

%% ILC parameters
gamma = 0.5; % 0.5, 0.85, 1.15, 1.5

%% AVAR inputs
p = floor(log2(number_of_iterations));
list_of_correlation_intervals = 2.^(0:p-2)'; % list of correlation intervals

%% Initialize variables to save MSE
tracking_error_mse = zeros(number_of_iterations,1);
true_error_mse     = zeros(number_of_iterations,1);

%% Initialize variables to save AVAR
vec_tracking_error = NaN(number_of_iterations,length_of_iteration);
vec_true_error     = NaN(number_of_iterations,length_of_iteration);
vec_tracking_avar  = NaN(p-1,length_of_iteration);
vec_true_avar      = NaN(p-1,length_of_iteration);

%% ILC Simulation
for index_run = 1:number_of_runs
    % Drift along iteration domain
    if 0<random_walk_coefficient
        random_walk = fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
            sampling_frequency,number_of_iterations); % Random walk
    else
        random_walk = zeros(number_of_iterations,1);
    end % NOTE: END IF statement
    
    for index_iter = 1:number_of_iterations
        %% Initialize variables
        plant_output      = NaN(length_of_iteration,1);
        true_plant_output = NaN(length_of_iteration,1);
        
        % Measurement noise along sample domain
        white_noise = fcn_AVAR_generateWhiteNoise(power_spectral_density,...
            sampling_frequency,length_of_iteration); % White noise
        
        % Control input update
        if 1==index_iter
            % Control input for first iteration is set equal to the reference input
            control_input = reference_input;
        else
            control_input = prev_control_input + ...
                gamma*[prev_tracking_error(2:length_of_iteration,1); ...
                       prev_tracking_error(length_of_iteration,1)];
        end % NOTE: END IF statement '1==index_iter'
        prev_control_input = control_input;
        
        for index_time = 1:length_of_iteration
            %% System in State-Space representation
            % State estimate: Discrete-time system
            if 1==index_time
                % Initial states are set such that the error at first sample is zero due to ICs
                system_state = (reference_input(1)/sum(C.^2))*C';
            else
                % State at '(index_time-1)*sampling_interval'
                system_state = A*system_state + B*control_input(index_time-1);
            end % NOTE: END IF statement '1==index_time'
            % Output measurement at '(index_time-1)*sampling_interval'
            plant_output(index_time,1)      = C*system_state + white_noise(index_time);
            true_plant_output(index_time,1) = C*system_state;
        end % NOTE: END FOR loop 'length_of_iteration'
        % Tracking error estimate
        delta_ref_input     = random_walk(index_iter)*ones(size(reference_input));
        tracking_error      = reference_input+delta_ref_input-plant_output;
        prev_tracking_error = tracking_error;
        true_tracking_error = reference_input-true_plant_output;
        tracking_error_mse(index_iter,1) = tracking_error_mse(index_iter,1) + mean(tracking_error.^2);
        true_error_mse(index_iter,1)     = true_error_mse(index_iter,1) + mean(true_tracking_error.^2);
        
        vec_tracking_error(index_iter,:) = tracking_error';
        vec_true_error(index_iter,:)     = true_tracking_error';
    end % NOTE: END FOR loop 'number_of_iterations'
end % NOTE: END FOR loop 'number_of_runs'
tracking_error_mse = tracking_error_mse/number_of_runs;
true_error_mse     = true_error_mse/number_of_runs;

%% Estimate AVAR in iteration domain
figure(12345)
clf
hold on
grid on
for index_time = 1:length_of_iteration
    vec_tracking_avar(:,index_time) = ...
        fcn_AVAR_favar([vec_tracking_error(:,index_time); 0], list_of_correlation_intervals);
    vec_true_avar(:,index_time) = ...
        fcn_AVAR_favar([vec_true_error(:,index_time); 0], list_of_correlation_intervals);
    plot(list_of_correlation_intervals, vec_true_avar(:,index_time))
end % NOTE: END FOR loop 'length_of_iteration'
set(gca,'YScale','log','XScale','log','FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Iterations]$','Interpreter','latex','FontSize',18)

% data.list_of_correlation_intervals = list_of_correlation_intervals;
% data.length_of_iteration = length_of_iteration;
% data.vec_tracking_avar = vec_tracking_avar;
% data.vec_true_avar     = vec_true_avar;
% save(file_name,"data");
% clear data
% load(file_name);
% 
% figure(12346)
% clf
% hold on
% grid on
% for index_time = 1:data.length_of_iteration
%     plot(data.list_of_correlation_intervals, data.vec_true_avar(:,index_time))
% end % NOTE: END FOR loop 'length_of_iteration'
% set(gca,'YScale','log','XScale','log','FontSize',13)
% ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
% xlabel('Correlation Interval $[Number \: of \: Iterations]$','Interpreter','latex','FontSize',18)
