clc; close all ; clear all;
load("data_set_mcv.mat")
load("data_set_sdlqr_wowind.mat")
load("straight_line1.mat")
data_set_lqr=data_set_sdlqr_wowind;

%% calculate varinace %%
traj_nom=traj_nom(1:11001);
n_ind=length(traj_nom); %points to be tracked
data_set=20; %num of sim


x_mcv=data_set_mcv(1:end,2);
y_mcv=data_set_mcv(1:end,3);
z_mcv=data_set_mcv(1:end,4);

x_mcv=reshape(x_mcv,[n_ind,data_set]);
y_mcv=reshape(y_mcv,[n_ind,data_set]);
z_mcv=reshape(z_mcv,[n_ind,data_set]);


x_lqr=data_set_lqr(1:data_set*n_ind,2);
y_lqr=data_set_lqr(1:data_set*n_ind,3);
z_lqr=data_set_lqr(1:data_set*n_ind,4);

%%
x_lqr=reshape(x_lqr,[n_ind,data_set]);
y_lqr=reshape(y_lqr,[n_ind,data_set]);
z_lqr=reshape(z_lqr,[n_ind,data_set]);

x_mcv_err=zeros(n_ind,data_set);y_mcv_err=x_mcv_err; z_mcv_err=x_mcv_err;
x_lqr_err=x_mcv_err;y_lqr_err=x_mcv_err;z_lqr_err=x_mcv_err;
for jj=1:data_set
    x_mcv_err(:,jj)= traj_nom(:,1)-x_mcv(:,jj);y_mcv_err(:,jj)= traj_nom(:,2)-y_mcv(:,jj);z_mcv_err(:,jj)= traj_nom(:,3)-z_mcv(:,jj);
    x_lqr_err(:,jj)= traj_nom(:,1)-x_lqr(:,jj);y_lqr_err(:,jj)= traj_nom(:,2)-y_lqr(:,jj);z_lqr_err(:,jj)= traj_nom(:,3)-z_lqr(:,jj);
    
end
%%

mean_traj_err_lqr=abs([ mean(x_lqr_err);mean(y_lqr_err);mean(z_lqr_err)]);
mean_traj_err_mcv=abs([ mean(x_mcv_err);mean(y_mcv_err);mean(z_mcv_err)]);
%%
% for jj=1:data_set
%     
%     x_rms_mcv(jj) = rms_error(x_mcv_err(:,jj)); x_var_mcv(jj)=var(x_mcv_err(:,jj));
%     y_rms_mcv(jj) = rms_error(y_mcv_err(:,jj)); y_var_mcv(jj)=var(y_mcv_err(:,jj));
%     z_rms_mcv(jj) = rms_error(z_mcv_err(:,jj)); z_var_mcv(jj)=var(z_mcv_err(:,jj));
%     
%     x_rms_lqr(jj) = rms_error(x_lqr_err(:,jj)); x_var_lqr(jj)=var(x_lqr_err(:,jj));
%     y_rms_lqr(jj) = rms_error(y_lqr_err(:,jj)); y_var_lqr(jj)=var(y_lqr_err(:,jj));
%     z_rms_lqr(jj) = rms_error(z_lqr_err(:,jj)); z_var_lqr(jj)=var(z_lqr_err(:,jj));
%     
% end

%%
for kk=1:n_ind
    
    x_v_mcv(kk)=var(x_mcv_err(kk,:)); x_v_mcv_max(kk)= max(x_mcv_err(kk,:));  x_v_mcv_min(kk)= min(x_mcv_err(kk,:)); 
    y_v_mcv(kk)=var(y_mcv_err(kk,:)); y_v_mcv_max(kk)= max(y_mcv_err(kk,:));  y_v_mcv_min(kk)= min(y_mcv_err(kk,:)); 
    z_v_mcv(kk)=var(z_mcv_err(kk,:)); z_v_mcv_max(kk)= max(z_mcv_err(kk,:));  z_v_mcv_min(kk)= min(z_mcv_err(kk,:)); 
    
    x_v_lqr(kk)=var(x_lqr_err(kk,:)); x_v_lqr_max(kk)= max(x_lqr_err(kk,:));  x_v_lqr_min(kk)= min(x_lqr_err(kk,:)); 
    y_v_lqr(kk)=var(y_lqr_err(kk,:));  y_v_lqr_max(kk)= max(y_lqr_err(kk,:));  y_v_lqr_min(kk)= min(y_lqr_err(kk,:)); 
    z_v_lqr(kk)=var(z_lqr_err(kk,:));  z_v_lqr_max(kk)= max(z_lqr_err(kk,:));  z_v_lqr_min(kk)= min(z_lqr_err(kk,:)); 
    
    x_rms_mcv(kk) = rms_error(x_mcv_err(kk,:));
    y_rms_mcv(kk) = rms_error(y_mcv_err(kk,:));
    z_rms_mcv(kk) = rms_error(z_mcv_err(kk,:));
    
    x_rms_lqr(kk) = rms_error(x_lqr_err(kk,:));
    y_rms_lqr(kk) = rms_error(y_lqr_err(kk,:));
    z_rms_lqr(kk) = rms_error(z_lqr_err(kk,:));
    
    
end

%%
%RMS ERROR %
sim=1:1:n_ind;
figure()
subplot(3,1,1)
plot(sim,x_rms_lqr,'--r','LineWidth',1.5)
hold on
plot(sim,x_rms_mcv,'-.b','LineWidth',2)
grid on

title("RMSE at each point over 20 simulation")
% xlabel("Time, [milisecond]")
ylabel("x-RMSE, [m]")
legend("LQR","MCV")
 xlim([0,11000])
%  xlim([0,1050])


subplot(3,1,2)
plot(sim,y_rms_lqr,'--r','LineWidth',1.5)
hold on
plot(sim,y_rms_mcv,'-.b','LineWidth',2)

grid on

% title("RMS error at each point over 50 simulation in y")
% xlabel("Time, [milisecond]")
ylabel("y-RMSE, [m]")
% legend("MCV","LQR")
 xlim([0,11000])
% xlim([0,1050])


subplot(3,1,3)
plot(sim,z_rms_lqr,'--r','LineWidth',1.5)
hold on
plot(sim,z_rms_mcv,'-.b','LineWidth',2)

grid on

% title("RMS error at each point over 50 simulation in z")
xlabel("Time, [milisecond]")
ylabel("z-RMSE, [m]")
 xlim([0,11000])
% legend("MCV","LQR")

%  xlim([0,1050])
%%
%Mean error %
sim=1:1:data_set;
figure()
plot(sim,mean_traj_err_mcv(1,:))
hold on
plot(sim,mean_traj_err_lqr(1,:))
title("Absolute mean trajectory error over 50 simulation in x")
xlabel("Simulation")
ylabel("Mean error, [m]")
legend("MCV","LQR")

figure()
plot(sim,mean_traj_err_mcv(2,:))
hold on
plot(sim,mean_traj_err_lqr(2,:))
title("Absolute mean trajectory error over 50 simulation in y")
xlabel("Simulation")
ylabel("Mean error, [m]")
legend("MCV","LQR")

figure()
plot(sim,mean_traj_err_mcv(3,:))
hold on
plot(sim,mean_traj_err_lqr(3,:))
title("Absolute mean trajectory error over 50 simulation in z")
xlabel("Simulation")
ylabel("Mean error, [m]")
legend("MCV","LQR")



%%
% Varinace %
n=1:1:n_ind;
figure()

subplot(3,1,1)
plot(n,x_v_lqr,'--r','LineWidth',1.5)
hold on
plot(n,x_v_mcv,'-.b','LineWidth',2)

grid on

title("Error variance {\sigma} at each point over 50 simulation")
% xlabel("Time, [milisecond]")
ylabel("x-{\sigma}, [m]")
legend("LQR","MCV")
xlim([0,1100])
%  xlim([0,1050])


subplot(3,1,2)
plot(n,y_v_lqr,'--r','LineWidth',1.5)
hold on
plot(n,y_v_mcv,'-.b','LineWidth',2)
grid on

% title("Error variance at each point over 50 simulation in y")
% xlabel("Time, [milisecond]")
ylabel("y-{\sigma}, [m]")
% legend("MCV","LQR")
 xlim([0,1100])
%  xlim([0,1050])

subplot(3,1,3)
plot(n,z_v_lqr,'--r','LineWidth',1.5)
hold on
plot(n,z_v_mcv,'-.b','LineWidth',2)
grid on

% title("Error variance at each point over 50 simulation in z")
xlabel("Time, [milisecond]")
ylabel("z-{\sigma}, [m]")
% legend("MCV","LQR")
 xlim([0,1100])
%  xlim([0,1050])
%%
figure()

id=length(traj_nom);
plot(data_set_lqr(1:id,1),data_set_lqr(1:id,2),'r','LineWidth',1.2)
hold on
% grid on
plot(data_set_mcv(1:id,1),data_set_mcv(1:id,2),'--b','LineWidth',1.2)
plot(data_set_mcv(1:id,1),traj_nom(:,1),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('x (m)')
title('Trajectory Tracking in LES wind')
legend("LQR","RSC","Ref")
% xlim([0,1400])


figure()
plot(data_set_lqr(1:id,1),data_set_lqr(1:id,3),'r','LineWidth',1.2)
hold on
% grid on
plot(data_set_mcv(1:id,1),data_set_mcv(1:id,3),'--b','LineWidth',1.2)
plot(data_set_mcv(1:id,1),traj_nom(:,2),'-.k','LineWidth',1.2)
% xlabel('time [s]')
ylabel('y (m)')
% title('trajectory in y direction')
% legend("LQR","Finite MCV","Nominal snap")
% xlim([0,1400])



figure()
plot(data_set_lqr(1:id,1),data_set_lqr(1:id,4),'r','LineWidth',1.2)
hold on
% grid on
plot(data_set_mcv(1:id,1),data_set_mcv(1:id,4),'--b','LineWidth',1.2)
plot(data_set_mcv(1:id,1),traj_nom(:,3),'-.k','LineWidth',1.2)

xlabel('times (s)')
ylabel('z (m)')
% xlim([0,1400])

% title('trajectory in z direction')
% legend("LQR","Finite MCV","Nominal snap")

%%
% figure()
% plot3(data_set_lqr(1:id,2),data_set_lqr(1:id,3),data_set_lqr(1:id,4),'--r','LineWidth',1.5)
% hold on 
% plot3(data_set_mcv(1:id,2),data_set_mcv(1:id,3),data_set_mcv(1:id,4),'-.b','LineWidth',2)
% plot3(traj_nom(1:id,1),traj_nom(1:id,2),traj_nom(1:id,3),'k')
% grid on
% xlabel("x [m]")
% ylabel("y [m]")
% zlabel("z [m]")
% legend("LQR","MCV","Nominal snap")

