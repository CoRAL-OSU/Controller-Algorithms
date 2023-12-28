clc
clear all
close all
%------ Main: State Dependent LQR without Wind information ------%
%                       Asma Tabassum                            %

%%
% Load required files
% Load nominal trajectory file, stable model gain K for the dynamic
% system ,variance matrices G and W
load("straight_line1.mat")  %calculated from minimum snap
load("W_mat.mat")
load("Gnw.mat")

%load("G_cov.mat")
% load("K_test.mat")
 load("LES_wind_updated.mat")


%%

% Simulation parameters
global uu kk
kk=0;
params.g = -9.81;
params.Q = diag([10,10,10,1,1,1,1,0.1,0.1,0.1]);
params.R = diag([1,5,5,0.1]);
params.G=G;
params.W=W;
params.gaama=0;  %this is highest
vw_bar=[2.72,1.752,0.0006];
v_zero=[0,0,0];
params.meanvel=v_zero;
%Constant wind
% wind =[0,0,-2];
% params.meanvel=wind;%[0;0;0];%WE DONT KNOW WIND
params.D=0*[0.218,0,0;
         0, 0.218 , 0;
         0 , 0, 0.004];
     
params.u0 = [0;0;0.02;7];  
params.x0 = [0;0;2;1;0;0;0;-params.meanvel.']';
params.xref=[0;0;2;1;0;0;0;-params.meanvel.']';

params.m=1.035 * 4/6; 
params.rho=1.15;
params.uref=[0;0;0;params.g*params.m];

%% find initial gain value
[A,B]=deriveLinSys_wowind(params.xref,params.u0,params);
K0=lqr(A,B,params.Q,params.R);
params.K=K0;
%% Data Preparation from Trajectory files


traj_pts=traj_nom(1:end,:);
traj_nom=traj_pts(:,2:11);
t_nom=traj_pts(:,1);
% u_nom=control(t_nom,traj_nom,params);



%% get the trajectory using new gains 

uu = [params.u0.'];
y0 = params.x0;

t0=0;

nomTraj=[t_nom, traj_nom];
final_time=t_nom(end);

Tf=1200; %total time
dt=0.01; 
X=params.x0;
params.wind=wind;
sim_num=1;
data_set=[];
for jj=1:sim_num
    1
%    wind_idx = (sim_num-1)*10+1;
%    wind_end=wind_idx+Tf/dt;
%    params.wind=wind(wind_idx:wind_end+1,:);
    [t,y]=sdef(Tf,dt,X,params,nomTraj);
    current_data=[t',y];
% 
%     data_set=[data_set;current_data];

end

%%
%saving datasets

    rpy=zeros(length(t),3);
    quat_angle=y(:,4:7);
for i=1:length(t)
    rpy(i,:)=quat2eul(quat_angle(i,:),"XYZ");
end

dataset=[t',y(:,1:3),y(:,8:10),rpy];

%% input data collect only for single run
u_set=uu(1:4:end,:);
gain_set=kk(1:4:end,:);
dataset=[dataset,u_set];

%% plot trajectory and input
traj_nom=params.xref.*ones(length(t),10);
figure()
plot(current_data(:,1),current_data(:,2))
hold on
plot(t,traj_nom(:,1))
xlabel('time [s]')
ylabel('x [m]')
title('trajectory in x direction')
legend('Actual','Nominal')

figure()
plot(current_data(:,1),current_data(:,3))
hold on
plot(t,traj_nom(:,2))
xlabel('time [s]')
ylabel('y [m]')
title('trajectory in y direction')
legend('Actual','Nominal')

figure()
plot(current_data(:,1),current_data(:,4))
hold on
plot(t,traj_nom(:,3))
xlabel('time [s]')
ylabel('z [m]')
title('trajectory in z direction')
legend('Actual','Nominal')

figure()
plot(current_data(:,1),u_set(:,1))
xlabel('time [s]')
ylabel('wx rad/s')
title('bodyrate in x direction, roll')

figure()
plot(current_data(:,1),u_set(:,2))
xlabel('time [s]')
ylabel('wy rad/s')
title('bodyrate in y direction, pitch')

figure()
plot(current_data(:,1),u_set(:,3))
xlabel('time [s]')
ylabel('wz rad/s')
title('bodyrate in z direction, yaw')

figure()
plot(current_data(:,1),u_set(:,4))
xlabel('time [s]')
ylabel('Cumulative Thrust, N')
title('Thrust')


%%
function u = finite_controller(y, params,iter)
global uu kk
p = reshape(y(1:3),[3 1]);
q = reshape(y(4:7),[4 1]);
v = reshape(y(8:10),[3 1]); %relative velocity
x=[p;q;v];


if iter==1
    prev_u=uu(end,:);%THIS SHOULD BE CONSTANT FOR ONE SIM so only taking the actual previous
else
    prev_u=uu(end-4,:);
end
% xdes=getref(t,nomTraj);
xdes=params.xref;
y;
[A,B] = deriveLinSys_wowind(y,prev_u,params);
K_now = lqr(A,B,params.Q,params.R);
kk=[kk;norm(K_now)];
x(8:10)=x(8:10)+params.meanvel.';%updating to track the ground velocity;
 u =params.u0 - K_now*(x-xdes');
 
 
 %limiting thrust
if u(4)>20
    u(4)=20;
end
upper_dist=[1;1;1;2.5];%0.4*u;%b
lower_dist=-[1;1;1;2.5];%-0.4*u;%a

udis=added_dist(lower_dist,upper_dist);
  u=udis+u;
end

%
%%
function udist= added_dist(a,b)

udist=a+(b-a).*rand(4,1);
end

%%
function [xdes] = getref(t,nomTraj)

% This function tracks time and match nominal trajectory at next nearest
% time
% Input:
% nomTraj is the nominal trajectory with t
% 

idxA = find(nomTraj(:,1)>t);
idx=idxA(1);
xdes=nomTraj(idx,2:end);

end


function dx = quad_dyn(y,params,nomTraj,i,iter)
% i is the time step coming from the sdef
global uu 
y;
p = reshape(y(1:3),[3 1]); %position
q = reshape(y(4:7),[4 1]);  %orinetation
v = reshape(y(8:10),[3 1]); % relative velocity
q=q/norm(q);
%LES wind
% wind=params.wind;
% id=uint32(i*100);
% wind_vel_next=wind(id+2,:);
% wind_vel_next(3)=-wind_vel_next(3);

%ground velocity to go into the controller
% if id==0 %at id=0; the y0 contain the ground velocity
%     vg = v;
% else
%     wind_vel=wind(id,:);
%     wind_vel(3)=-wind_vel(3);%positive up data
%     vg= v+wind_vel';
% end


%Constant wind

% wind_vel=params.meanvel;
% vg= v+wind_vel';%vg=vr+vw;
% wind_vel_next=wind_vel;

%LES wind
% wind=params.wind;
% id=uint32(i*100);
% wind_vel_next=wind(id+1,:);


u = finite_controller(y, params,iter);
uu=[uu;u'];

wn = reshape(u(1:3),[3 1]);
cn = [0;0;u(4)];

% G=  params.G;
% W=  params.W;
% noise=G*chol(W)*G'*rand(10,1);

dp=double(v);%+wind_vel_next';
dq = double(0.5*quater_mult(q)*[0;wn]);
dvc = double([0;0;params.m*params.g] + quater_rot(q)*cn);
R=quater_rot(q);
dvd=  0.5*params.rho*R*params.D*R'*norm(-v)*-v;
dv=(dvc+dvd)/params.m;


dx = [dp;dq;dv];


dx=dx';


end


function R = quater_rot(q)

w = q(1);
x = q(2);
y = q(3);
z = q(4);

Q = [w -x -y -z;...
     x w -z y;...
     y z w -x;...
     z -y x w];

Qc =[w -x -y -z;...
     x w z -y;...
     y -z w x;...
     z y -x w];
t = Qc'*Q;
R = t(2:4,2:4);
end


function Q = quater_mult(q)

w = q(1);
x = q(2);
y = q(3);
z = q(4);

Q = [w -x -y -z;...
     x w -z y;...
     y z w -x;...
     z -y x w];

end

%% Helper functions


function [t,X]=sdef(Tf,dt,X,params,nomTraj)  % Runge Kutta 4
t = 0:dt:Tf;
iter=0;
for i=1:length(t)-1    
    iter=+1;
    X(i,:);
    f1 = quad_dyn(X(i,:),params,nomTraj,t(i),iter) ;
    f2 = quad_dyn(X(i,:)+dt/2*f1,params,nomTraj,t(i),iter);
    f3 = quad_dyn(X(i,:)+dt/2*f2,params,nomTraj,t(i),iter);
    f4 = quad_dyn(X(i,:)+dt*f3,params,nomTraj,t(i),iter);
    X(i+1,:)=X(i,:)+(dt/6)*(f1+2*f2+2*f3+f4);
end
end

