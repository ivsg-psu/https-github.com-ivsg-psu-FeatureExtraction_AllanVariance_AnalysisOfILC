%%%%%%%%%%%%%%%%%%% script_plot_InDrift_OutWhiteNoise.m %%%%%%%%%%%%%%%%%%%
%% Purpose:
%   To plot MSE (Both True and Tracking) for demonstrating the effect of 
%   1,2--> Drift (Random Walk) on the performance of ILC tracking a drifting PRBS
%   3,4--> Measurement noise (White) on the performance of ILC tracking a drifting PRBS
%
% Author:  Satya Prasad
% Updated: 2025/02/10

%% Prepare the workspace
clearvars
close all
clc

%% Add path to dependencies
addpath('../Data')

%% Plot Drift variation
figure(01)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
load("20250209_InDrift0p0_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','Color','#FF8800','LineWidth',1.2)
load("20250209_InDrift0p01_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','b','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','r','LineWidth',1.2)
load("20250209_InDrift0p1_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','g','LineWidth',1.2)
load("20250209_InDrift0p2_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','k','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
legend('No Drift','$C_{rw} = 0.01$','$C_{rw} = 0.02$','$C_{rw} = 0.1$','$C_{rw} = 0.2$',...
    'Interpreter','latex','Location','best','FontSize',13)

figure(02)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
load("20250209_InDrift0p0_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','Color','#FF8800','LineWidth',1.2)
load("20250209_InDrift0p01_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','b','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','r','LineWidth',1.2)
load("20250209_InDrift0p1_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','g','LineWidth',1.2)
load("20250209_InDrift0p2_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','k','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
set(gca,'XScale','log','XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
legend('No Drift','$C_{rw} = 0.01$','$C_{rw} = 0.02$','$C_{rw} = 0.1$','$C_{rw} = 0.2$',...
    'Interpreter','latex','Location','best','FontSize',13)
ylim([0.02 0.1])

%% Plot Noise variation
figure(03)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','r','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p00004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','g','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p000004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','k','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
legend('$C_{wn} = 0.0004$','$C_{wn} = 0.00004$','$C_{wn} = 0.000004$',...
    'Interpreter','latex','Location','best','FontSize',13)

figure(04)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','r','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p00004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','g','LineWidth',1.2)
load("20250209_InDrift0p02_OutWhite0p000004.mat")
plot(1:data.number_of_iterations,data.tracking_error_mse','k','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
legend('$C_{wn} = 0.0004$','$C_{wn} = 0.00004$','$C_{wn} = 0.000004$',...
    'Interpreter','latex','Location','best','FontSize',13)

%% Reset vs Classical ILC
figure(11)
clf
width = 540; height = 400; right = 100; bottom = 100;
set(gcf, 'position', [right, bottom, width, height])
hold on
grid on
load("20250209_InDrift0p02_OutWhite0p0004.mat")
plot(1:data.number_of_iterations,data.true_error_mse','r','LineWidth',2.4)
load("20250218_InDrift0p02_OutWhite0p0004_Reset.mat")
plot(1:data.number_of_iterations,data.true_error_mse','g','LineWidth',1.2)
set(gca,'YScale','log','XScale','log','XTick',[1e0 1e1 1e2 1e3 1e4 1e5],'FontSize',13)
ylabel('Mean Squared Error $[Unit^{2}]$','Interpreter','latex','FontSize',18)
xlabel('Iteration Number','Interpreter','latex','FontSize',18)
legend('Classical ILC','Reset ILC','Interpreter','latex','Location','best','FontSize',13)
ylim([1e-3 1e0])
