clc; close all ; clear all;
%compare tracking and input for state dependent LQR 

load("stress_testing1_lqrwind_n.mat")
load("stress_testing1_lqrwowind_n.mat")
load("straight_line1.mat")
t=traj_nom(1:11001,1);
traj_nom=traj_nom(1:11001,2:11);

drag=stress_testing1_lqrwind_n;
no_drag=stress_testing1_lqrwowind_n;

%%
id=11001;
figure()

subplot(3,1,1)

plot(t,drag(:,2),'r','LineWidth',1.2)
hold on
grid on
plot(t,no_drag(:,2),'--b','LineWidth',1.2)
plot(t,traj_nom(:,1),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('x (m)')
% title('Trajectory Tracking in constant wind [5,4,1]m')
legend(" SD LQR with Drag model","SD LQR without Drag model","Nominal snap")
% xlim([0,1400])
title('Trajectory Tracking in constant wind [5,4,1]m')


subplot(3,1,2)
plot(t,drag(:,3),'r','LineWidth',1.2)
hold on
grid on
plot(t,no_drag(:,3),'--b','LineWidth',1.2)
plot(t,traj_nom(:,2),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('y (m)')
% title('trajectory in y direction')
% legend(" SD LQR with Drag model","SD LQR without Drag model","Nominal snap")
% xlim([0,1400])
% title('Trajectory Tracking in constant wind [5,4,1]m')


subplot(3,1,3)
plot(t,drag(:,4),'r','LineWidth',1.2)
hold on
 grid on
plot(t,no_drag(:,4),'--b','LineWidth',1.2)
plot(t,traj_nom(:,3),'-.k','LineWidth',1.2)
% legend(" SD LQR with Drag model","SD LQR without Drag model","Nominal snap")
xlabel('times (s)')
ylabel('z (m)')
% title('Trajectory Tracking in constant wind [5,4,1]m')

%% Input

figure()
subplot(2,2,1)


plot(t,drag(:,12),'r','LineWidth',1.2)
hold on
grid on
plot(t,no_drag(:,12),'--b','LineWidth',1.2)

% xlabel('time [s]')
ylabel('wx rad/s')
legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
% xlim([0,1400])

title('Input in constant wind [5,4,1]m')

subplot(2,2,2)
plot(t,drag(:,13),'r','LineWidth',1.2)
hold on
grid on
plot(t,no_drag(:,13),'--b','LineWidth',1.2)
% xlabel('time [s]')
ylabel('wy rad/s')
% title('trajectory in y direction')
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
% % xlim([0,1400])
% title('Input in constant wind [5,4,1]m')


subplot(2,2,3)
plot(t,drag(:,14),'r','LineWidth',1.2)
hold on
grid on
plot(t,no_drag(:,14),'--b','LineWidth',1.2)

% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
xlabel('times (s)')
ylabel('wz rad/s')

% title('Input in constant wind [5,4,1]m')


subplot(2,2,4)
plot(t,drag(:,15),'r','LineWidth',1.2)
hold on
 grid on
plot(t,no_drag(:,15),'--b','LineWidth',1.2)

% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
xlabel('times (s)')
ylabel('fc N')

% title('Input in constant wind [5,4,1]m')