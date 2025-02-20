%%%%%%%%%%%%%%% script_classicalILC_InDrift_OutWhiteNoise.m %%%%%%%%%%%%%%%
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
number_of_iterations = 2^15;
number_of_runs = 100;
flag_input = 2; % 1 -> Sinusoid input, 2 -> PRBS input
reset_iteration = 128;

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

%% Initialize variables to save MSE
tracking_error_mse = zeros(number_of_iterations,1);
true_error_mse     = zeros(number_of_iterations,1);

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
        if mod(index_iter,reset_iteration)==1 % Resetting drift after every 'reset_iteration' iterations
            random_walk(index_iter:number_of_iterations) = ...
                fcn_AVAR_generateRandomWalk(random_walk_coefficient,...
                sampling_frequency,number_of_iterations-index_iter+1); % Random walk
        end % NOTE: END IF statement
        
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
        elseif mod(index_iter,reset_iteration)==1
            % Control input is same as that of the previous iteration
            control_input = reference_input;
%             continue
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
    end % NOTE: END FOR loop 'number_of_iterations'
end % NOTE: END FOR loop 'number_of_runs'
tracking_error_mse = tracking_error_mse/number_of_runs;
true_error_mse     = true_error_mse/number_of_runs;

%% Plot the results
figure(01)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(1:number_of_iterations,true_error_mse','b','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)

figure(02)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
plot(1:number_of_iterations,tracking_error_mse','b','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
