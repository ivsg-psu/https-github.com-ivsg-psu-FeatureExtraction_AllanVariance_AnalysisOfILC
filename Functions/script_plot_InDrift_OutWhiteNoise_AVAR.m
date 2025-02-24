%%%%%%%%%%%%%%%%% script_plot_InDrift_OutWhiteNoise_AVAR.m %%%%%%%%%%%%%%%%
%% Purpose:
%   To plot MSE (Both True and Tracking) for demonstrating the effect of 
%   1--> Drift (Random Walk) on resetting ILC
%   2--> Measurement noise (White) on resetting ILC
%
% Author:  Satya Prasad
% Updated: 2025/02/18

%% Prepare the workspace
clearvars
close all
clc

%% Add path to dependencies
addpath('../Data')

%% Define inputs and other parameters
power_spectral_density  = 0.0004;

%% Plot Drift variation
figure(01)
clf
width = 540; height = 770.04; right = 100; bottom = 10;
set(gcf, 'position', [right, bottom, width, height])
axis_position = [0.1354, 432.72/height, 0.7696, 307.32/height];
subplot(2,1,1)
hold on
grid on
load("20250209_InDrift0p01_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','b','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','r','LineWidth',1.2)
load("20250209_InDrift0p1_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','g','LineWidth',1.2)
load("20250209_InDrift0p2_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','k','LineWidth',1.2)
legend('$C_{rw} = 0.01$','$C_{rw} = 0.02$','$C_{rw} = 0.1$','$C_{rw} = 0.2$',...
    'Interpreter','latex','Location','best','FontSize',13)
set(gca,'Position',axis_position,'YScale','log','XScale','log',...
    'XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([1e-3 1e2])

axis_position = [0.1354, 62.7/height, 0.7696, 307.32/height];
subplot(2,1,2)
hold on
grid on
load("20250210_InDrift0p01_OutWhite0p0004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.01,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'b')
end % NOTE: END FOR loop 'length_of_iteration'
load("20250210_InDrift0p02_OutWhite0p0004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.02,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'r')
end % NOTE: END FOR loop 'length_of_iteration'
load("20250210_InDrift0p1_OutWhite0p0004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.1,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'g')
end % NOTE: END FOR loop 'length_of_iteration'
load("20250210_InDrift0p2_OutWhite0p0004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.2,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'k')
end % NOTE: END FOR loop 'length_of_iteration'
set(gca,'Position',axis_position,'YScale','log','XScale','log',...
    'XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Iterations]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
ylim([10^-4.5 1e1])

%% Plot Noise variation
figure(02)
clf
width = 540; height = 770.04; right = 100; bottom = 10;
set(gcf, 'position', [right, bottom, width, height])
axis_position = [0.1354, 432.72/height, 0.7696, 307.32/height];
subplot(2,1,1)
hold on
grid on
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','r','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p00004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','g','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p000004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','k','LineWidth',1.2)
legend('$C_{wn} = 0.0004$','$C_{wn} = 0.00004$','$C_{wn} = 0.000004$',...
    'Interpreter','latex','Location','best','FontSize',13)
set(gca,'Position',axis_position,'YScale','log','XScale','log',...
    'XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
title('$(a)$','Interpreter','latex','FontSize',18)
ylim([1e-4 1e0])

axis_position = [0.1354, 62.7/height, 0.7696, 307.32/height];
subplot(2,1,2)
hold on
grid on
load("20250210_InDrift0p02_OutWhite0p0004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.02,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'r')
end % NOTE: END FOR loop 'length_of_iteration'
load("20250210_InDrift0p02_OutWhite0p00004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.02,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'g')
end % NOTE: END FOR loop 'length_of_iteration'
load("20250210_InDrift0p02_OutWhite0p000004_AVAR.mat")
avar_calculated = fcn_AVAR_avarRandomWalk(0.02,data.list_of_correlation_intervals/50);
for index_time = 2:data.length_of_iteration
    plot(data.list_of_correlation_intervals,...
        data.vec_tracking_avar(:,index_time)+avar_calculated,'k')
end % NOTE: END FOR loop 'length_of_iteration'
set(gca,'Position',axis_position,'YScale','log','XScale','log',...
    'XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Allan Variance $[Unit^2]$','Interpreter','latex','FontSize',18)
xlabel('Correlation Interval $[Number \: of \: Iterations]$','Interpreter','latex','FontSize',18)
title('$(b)$','Interpreter','latex','FontSize',18)
