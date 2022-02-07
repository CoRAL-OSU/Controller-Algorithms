clc; close all ; clear all;
%compare tracking and input for state dependent LQR 

load("les_lqrwind_1")
load("les_lqrwind_2")
load("les_lqrwind_3")
load("les_lqrwowind_1")
load("les_lqrwowind_2")
load("les_lqrwowind_3")
load("straight_line1.mat")
%%
t=traj_nom(1:11001,1);
traj_nom=traj_nom(1:11001,2:11);

drag=les_lqrwind_1;
no_drag=les_lqrwowind_1;
drag2=les_lqrwind_2;
no_drag2=les_lqrwowind_3;
drag3=les_lqrwind_3;
no_drag3=les_lqrwowind_3;

%%
id=11001;
figure()

subplot(3,1,1)
% 
% plot(t,drag(:,2),'r','LineWidth',1.2)
% hold on
% plot(t,no_drag(:,2),'--b','LineWidth',1.2)
plot(t,drag2(:,2),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,2),'b','LineWidth',1.2)
plot(t,drag3(:,2),'r:','LineWidth',1.2)
plot(t,no_drag3(:,2),'b:','LineWidth',1.2)
plot(t,traj_nom(:,1),'-.k','LineWidth',1.2)
ylabel('x (m)')
grid on
% title('Trajectory Tracking in constant wind [5,4,1]m')
% xlim([0,1400])
title('Trajectory Tracking in LES')


subplot(3,1,2)
% plot(t,drag(:,3),'r','LineWidth',1.2)
% hold on
% plot(t,no_drag(:,3),'--b','LineWidth',1.2)

plot(t,drag2(:,3),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,3),'b','LineWidth',1.2)
plot(t,drag3(:,3),'r:','LineWidth',1.2)
plot(t,no_drag3(:,3),'b:','LineWidth',1.2)
plot(t,traj_nom(:,2),'-.k','LineWidth',1.2)
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
plot(t,drag2(:,4),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,4),'b','LineWidth',1.2)
plot(t,drag3(:,4),'r:','LineWidth',1.2)
plot(t,no_drag3(:,4),'b:','LineWidth',1.2)
plot(t,traj_nom(:,3),'-.k','LineWidth',1.2)
% legend(" SD LQR with Drag model","SD LQR without Drag model","Nominal snap")
xlabel('times (s)')
grid on
ylabel('z (m)')
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Model 3:SD LQR with Drag model","Model 3: SD LQR without Drag model","Nominal")
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Nominal")
legend("Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Model 3:SD LQR with Drag model","Model 3: SD LQR without Drag model","Nominal")

% title('Trajectory Tracking in constant wind [5,4,1]m')

%% Input

figure()
subplot(2,2,1)


% plot(t,drag(:,12),'r','LineWidth',1.2)
% hold on
% grid on
% plot(t,no_drag(:,12),'--b','LineWidth',1.2)
plot(t,drag2(:,12),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,12),'b','LineWidth',1.2)

plot(t,drag3(:,12),'r:','LineWidth',1.2)
plot(t,no_drag3(:,12),'b:','LineWidth',1.2)
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
plot(t,drag2(:,13),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,13),'b','LineWidth',1.2)
plot(t,drag3(:,13),'r:','LineWidth',1.2)
plot(t,no_drag3(:,13),'b:','LineWidth',1.2)
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
plot(t,drag2(:,14),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,14),'b','LineWidth',1.2)
plot(t,drag3(:,14),'r:','LineWidth',1.2)
plot(t,no_drag3(:,14),'b:','LineWidth',1.2)
grid on
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
xlabel('times (s)')
ylabel('wz rad/s')
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Model 3:SD LQR with Drag model","Model 3: SD LQR without Drag model")
% legend("Model 1:SD LQR with Drag model","Model 1: SD LQR without Drag model","Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Nominal")
legend("Model 2:SD LQR with Drag model","Model 2: SD LQR without Drag model","Model 3:SD LQR with Drag model","Model 3: SD LQR without Drag model","Nominal")


% title('Input in constant wind [5,4,1]m')


subplot(2,2,4)
% plot(t,drag(:,15),'r','LineWidth',1.2)
% hold on
%  grid on
% plot(t,no_drag(:,15),'--b','LineWidth',1.2)
plot(t,drag2(:,15),'--r','LineWidth',1.2)
hold on
plot(t,no_drag2(:,15),'b','LineWidth',1.2)
plot(t,drag3(:,15),'r:','LineWidth',1.2)
plot(t,no_drag3(:,15),'b:','LineWidth',1.2)
grid on
% legend(" SD LQR with Drag model","SD LQR without Drag model")%,"Nominal snap")
xlabel('times (s)')
ylabel('fc N')

% title('Input in constant wind [5,4,1]m')