%%%%%%%%%%%%%% script_classicalILC_InZeroDrift_OutZeroNoise.m %%%%%%%%%%%%%
%% Purpose:
%   To demonstrate ILC using the example from 
%   K. L. Moore, Y. Chen and H. -S. Ahn, "Iterative Learning Control: 
%   A Tutorial and Big Picture View," Proceedings of the 45th IEEE 
%   Conference on Decision and Control, San Diego, CA, USA, 2006.
%
%   Input: Sinusoid or PRBS (Pseudo Random Binary Sequence)
%   Measurement Noise: None
%
% Author:  Satya Prasad
% Updated: 2025/02/10

%% Prepare the workspace
clearvars
close all
clc

%% Define inputs and other parameters
number_of_iterations = 20;
flag_input = 2; % 1 -> Sinusoid input, 2 -> PRBS input

sampling_frequency  = 50; % Sampling frequency [Hz]
sampling_interval   = 1/sampling_frequency; % Sampling interval/Simulation step size [s]
start_time          = 0; % Seconds
end_time            = 100*(2*pi*(100/8)*sampling_interval); % Seconds
time_vector         = start_time:sampling_interval:end_time;
length_of_iteration = length(time_vector); % Length for which one iteration lasts
dis_time_vector     = 0:length_of_iteration-1;
iterations_to_plot  = [1, 2, 4, 7, 20];

%% Plant definition
A = [-0.8, -0.22;...
        1, 0];
B = [0.5, 1]';
C = [1, 0.5];
D = 0;

%% Desired trajectory
if 1==flag_input
    reference_input = sin(8/100*dis_time_vector); % Desired trajectory to be tracked
elseif 2==flag_input
    rng('default')
    coin_toss   = rand(length_of_iteration,1)';
    flip_sign   = coin_toss>.95;
    total_flips = cumsum(flip_sign);
    reference_input = mod(total_flips,2); % Desired trajectory to be tracked
end % NOTE: END IF statement 'flag_input'

%% ILC parameters
list_of_gamma   = [0.5, 0.85, 1.15, 1.5]; % Tunable parameter for ILC error (Learning gain)
number_of_gamma = numel(list_of_gamma);

%% Initialize variables
tracking_error = NaN(number_of_iterations,length_of_iteration,number_of_gamma);
control_input  = NaN(number_of_iterations,length_of_iteration,number_of_gamma);
plant_output   = NaN(number_of_iterations,length_of_iteration,number_of_gamma);
error_mse      = NaN(number_of_iterations,number_of_gamma);

%% ILC Simulation
for index_gamma = 1:number_of_gamma
    gamma = list_of_gamma(index_gamma);
    
    for index_iter = 1:number_of_iterations
        % Control input update
        if 1==index_iter
            % Control input for first iteration is set equal to the reference input
            control_input(index_iter,:,index_gamma) = reference_input;
        else
            control_input(index_iter,:,index_gamma) = ...
                control_input(index_iter-1,:,index_gamma) + ...
                gamma*[tracking_error(index_iter-1,2:length_of_iteration,index_gamma), ...
                       tracking_error(index_iter-1,length_of_iteration,index_gamma)];
        end % NOTE: END IF statement '1==index_iter'
        
        for index_time = 1:length_of_iteration
            %% System in State-Space representation
            % State estimate: Discrete-time system
            if 1==index_time
                % Initial states are set such that the error at first sample is zero due to ICs
                system_state = (reference_input(1)/sum(C.^2))*C';
            else
                % State at '(index_time-1)*sampling_interval'
                system_state = A*system_state + B*control_input(index_iter,index_time-1,index_gamma);
            end % NOTE: END IF statement '1==index_time'
            % Output measurement
            plant_output(index_iter,index_time,index_gamma) = C*system_state; % Output at '(index_time-1)*sampling_interval'
        end % NOTE: END FOR loop 'length_of_iteration'
        % Tracking error estimate
        tracking_error(index_iter,:,index_gamma) = reference_input-plant_output(index_iter,:,index_gamma);
        error_mse(index_iter,index_gamma)        = mean(tracking_error(index_iter,:,index_gamma).^2);
    end % NOTE: END FOR loop 'number_of_iterations'
end % NOTE: END FOR loop 'number_of_gamma'

%% Plot the results
number_of_plots  = numel(iterations_to_plot);
custom_color_map = {'#FF8800','b','k','g','r'};
%%% Plant output: One plot per learning gain
for index_gamma = 1:number_of_gamma
    legend_cell = cell(number_of_plots+1,1);
    figure(index_gamma)
    clf
    width = 540; height = 400; right = 100; bottom = 100;
    set(gcf, 'position', [right, bottom, width, height])
    hold on
    grid on
    plot(dis_time_vector,reference_input,'Color',[0.8, 0.8, 0.8],'LineWidth',5)
    legend_cell{1} = "Desired";
    for index_plot = 1:number_of_plots
        index_iter = iterations_to_plot(index_plot);
        plot(dis_time_vector,plant_output(index_iter,:,index_gamma),...
            'Color',custom_color_map{index_plot},'LineWidth',1.2)
        legend_cell{index_plot+1} = strcat("iter ", num2str(index_iter));
    end % NOTE: END FOR loop 'number_of_plots'
    legend(legend_cell,'NumColumns',3,'Location','best','Interpreter','latex','FontSize',13)
    set(gca,'FontSize',13)
    ylabel('Plant Output $[Unit]$','Interpreter','latex','FontSize',18)
    xlabel('Time Step $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
    title(['$\gamma = $ ' num2str(list_of_gamma(index_gamma))],'Interpreter','latex','FontSize',18)
    xlim([0 200])
    ylim([-0.5 1.5])
end % NOTE: END FOR loop 'number_of_gamma'

%%% Tracking error: One plot per learning gain
for index_gamma = 1:number_of_gamma
    legend_cell = cell(number_of_plots,1);
    figure(number_of_gamma+index_gamma)
    clf
    width = 540; height = 400; right = 100; bottom = 100;
    set(gcf, 'position', [right, bottom, width, height])
    hold on
    grid on
    for index_plot = 1:number_of_plots
        index_iter = iterations_to_plot(index_plot);
        plot(dis_time_vector,tracking_error(index_iter,:,index_gamma),...
            'Color',custom_color_map{index_plot},'LineWidth',1.2)
        legend_cell{index_plot} = strcat("iter ", num2str(index_iter));
    end % NOTE: END FOR loop 'number_of_plots'
    legend(legend_cell,'NumColumns',3,'Location','best','Interpreter','latex','FontSize',13)
    set(gca,'FontSize',13)
    ylabel('Tracking Error $[Unit]$','Interpreter','latex','FontSize',18)
    xlabel('Time Step $[Number \: of \: Samples]$','Interpreter','latex','FontSize',18)
    title(['$\gamma = $ ' num2str(list_of_gamma(index_gamma))],'Interpreter','latex','FontSize',18)
    xlim([0 200])
end % NOTE: END FOR loop 'number_of_gamma'

%%% MSE
custom_color_map = {'r','#FF8800','b','g'};
legend_cell = cell(number_of_gamma,1);
figure(2*number_of_gamma+1)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
for index_gamma = 1:number_of_gamma
    plot(1:number_of_iterations,error_mse(:,index_gamma)',...
        'Color',custom_color_map{index_gamma},'LineWidth',1.2)
    legend_cell{index_gamma} = ['$\gamma = $ ' num2str(list_of_gamma(index_gamma))];
end % NOTE: END FOR loop 'number_of_gamma'
legend(legend_cell,'NumColumns',2,'Location','best','Interpreter','latex','FontSize',13)
set(gca,'YScale','log','FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
