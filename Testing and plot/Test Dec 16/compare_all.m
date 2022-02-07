clc; close all ; clear all;
%compare tracking and input for state dependent LQR 

load("les_lqrwind_1.mat")
load("les_offlqrwind.mat")
load("les_offmcvwind.mat")


load("straight_line1.mat")
%%
t=traj_nom(1:11001,1);
traj_nom=traj_nom(1:11001,2:11);


sdlqr=les_lqrwind_1;
offlqr=les_offlqrwind;
offmcv=les_offmcvwind;


%%
id=20/0.01;
t=t(1:1:id);
figure()

subplot(3,1,1)
plot(t,sdlqr(1:id,2),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,2),'r:','LineWidth',1.2)
plot(t,offmcv(1:id,2),'b','LineWidth',1.2)
plot(t,traj_nom(1:id,1),'-.k','LineWidth',1.2)
ylabel('x (m)')
grid on
% title('Trajectory Tracking in constant wind [5,4,1]m')
% xlim([0,1400])
title('Trajectory Tracking in LES')


subplot(3,1,2)
% plot(t,drag(:,3),'r','LineWidth',1.2)
% hold on
% plot(t,no_drag(:,3),'--b','LineWidth',1.2)

plot(t,sdlqr(1:id,3),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,3),'r:','LineWidth',1.2)
plot(t,offmcv(1:id,3),'b','LineWidth',1.2)
plot(t,traj_nom(1:id,2),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('y (m)')
grid on
% title('trajectory in y direction')
% legend(" SD LQR with Drag model","SD LQR without Drag model","Nominal snap")
% xlim([0,1400])
% title('Trajectory Tracking in constant wind [5,4,1]m')


subplot(3,1,3)
% plot(t,drag(:,4),'r','LineWidth',1.2)
% hold on
%  grid on
% plot(t,no_drag(:,4),'--b','LineWidth',1.2)
plot(t,sdlqr(1:id,4),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,4),'r:','LineWidth',1.2)
plot(t,offmcv(1:id,4),'b','LineWidth',1.2)

plot(t,traj_nom(1:id,3),'-.k','LineWidth',1.2)
% legend(" SD LQR with Drag model","SD LQR without Drag model","Nominal snap")
xlabel('times (s)')
grid on
ylabel('z (m)')
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Model 3:SD LQR with Drag model","Model 3: SD LQR without Drag model","Nominal")
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Nominal")
legend("SD LQR","Offline LQR","Offline MCV","Nominal")

% title('Trajectory Tracking in constant wind [5,4,1]m')

%% Input

figure()
subplot(2,2,1)


% plot(t,drag(:,12),'r','LineWidth',1.2)
% hold on
% grid on
% plot(t,no_drag(:,12),'--b','LineWidth',1.2)
plot(t,sdlqr(1:id,12),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,12),'r:','LineWidth',1.2)

plot(t,offmcv(1:id,12),'b','LineWidth',1.2)
ylabel('wx rad/s')
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
% xlim([0,1400])
grid on
title('Input in LES')

subplot(2,2,2)
% plot(t,drag(:,13),'r','LineWidth',1.2)
% hold on
% grid on
% plot(t,no_drag(:,13),'--b','LineWidth',1.2)
plot(t,sdlqr(1:id,13),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,13),'r:','LineWidth',1.2)
plot(t,offmcv(1:id,13),'b','LineWidth',1.2)
% xlabel('time [s]')
ylabel('wy rad/s')
grid on
% title('trajectory in y direction')
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
% % xlim([0,1400])
% title('Input in constant wind [5,4,1]m')


subplot(2,2,3)
% plot(t,drag(:,14),'r','LineWidth',1.2)
% hold on
% grid on
% plot(t,no_drag(:,14),'--b','LineWidth',1.2)
plot(t,sdlqr(1:id,14),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,14),'r:','LineWidth',1.2)
plot(t,offmcv(1:id,14),'b','LineWidth',1.2)
grid on
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
xlabel('times (s)')
ylabel('wz rad/s')
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Model 3:SD LQR with Drag model","Model 3: SD LQR without Drag model")
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Nominal")
legend("SD LQR","Offline LQR","Offline MCV")

% title('Input in constant wind [5,4,1]m')


subplot(2,2,4)
% plot(t,drag(:,15),'r','LineWidth',1.2)
% hold on
%  grid on
% plot(t,no_drag(:,15),'--b','LineWidth',1.2)
plot(t,sdlqr(1:id,15),'--r','LineWidth',1.2)
hold on
plot(t,offlqr(1:id,15),'r:','LineWidth',1.2)
plot(t,offmcv(1:id,15),'b','LineWidth',1.2)

grid on
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
xlabel('times (s)')
ylabel('fc N')

% title('Input in constant wind [5,4,1]m')