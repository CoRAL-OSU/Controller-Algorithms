clc; close all ; clear all;
load("data_set_mcv.mat")
load("data_set_sdlqr_wowind.mat")
load("data_set_lqr.mat")
load("data_set_sdlqr.mat")
load("straight_line1.mat")
traj_nom=traj_nom(1:11001);

load("data_set_simulink_sdlqr.mat")
load("data_set_simulink_mcv.mat")
%%

id=11001;
plot(data_set_lqr(1:id,1),data_set_sdlqr_wowind(1:id,2),'r','LineWidth',1.2)
hold on
grid on
plot(data_set_mcv(1:id,1),data_set_mcv(1:id,2),'--b','LineWidth',1.2)
plot(data_set_mcv(1:id,1),traj_nom(:,1),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('x (m)')
title('Trajectory Tracking in LES wind')
legend(" SD LQR WO WIND","Finite MCV","Nominal snap")
% xlim([0,1400])


figure()
plot(data_set_lqr(1:id,1),data_set_sdlqr_wowind(1:id,3),'r','LineWidth',1.2)
hold on
% grid on
plot(data_set_mcv(1:id,1),data_set_mcv(1:id,3),'--b','LineWidth',1.2)
plot(data_set_mcv(1:id,1),traj_nom(:,2),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('y (m)')
% title('trajectory in y direction')
legend(" SD LQR WO WIND","Finite MCV","Nominal snap")
% xlim([0,1400])



figure()
plot(data_set_lqr(1:id,1),data_set_sdlqr_wowind(1:id,4),'r','LineWidth',1.2)
hold on
% grid on
plot(data_set_mcv(1:id,1),data_set_mcv(1:id,4),'--b','LineWidth',1.2)
plot(data_set_mcv(1:id,1),traj_nom(:,3),'-.k','LineWidth',1.2)
legend(" SD LQR WO WIND","Finite MCV","Nominal snap")
xlabel('times (s)')
ylabel('z (m)')