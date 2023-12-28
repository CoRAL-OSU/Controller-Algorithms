clc
clear all
close all
%------ Main: State Dependent LQR with Wind information ------%
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
vw_bar=[2.72,1.752,0.0006];
params.meanvel=vw_bar;
%Constant wind
% wind =[0,0,-2];
%params.meanvel=wind;


params.D=[0.318,0,0;
         0, 0.318 , 0;
         0 , 0, 0.004];
params.m=1.035 * 4/6; 
params.rho=1.15; 
params.u0 = [0;0;0.02;7]; 
params.x0 = [0;0;2;1;0;0;0;-params.meanvel.';0;0;0]';
params.xref_l=[0;0;2;1;0;0;0;-params.meanvel.']';
params.l=0.225 ;                                           % arm length [m]
params.Ix= 0.0469 ;
params.Iy= 0.0358 ;
params.Iz=  0.101 * 4/6;
params.b =  1.5652e-08;% * (60/(2*pi))^2 ;                  % Thrust coeffcient intially N/(rpm^2), %now N/(rev/s^2)
params.k =  2.0862e-10;%* (60/(2*pi))^2 ;   

%% find initial gain value
[A,B]=deriveLinSys(params.xref_l,params.u0,params);
K0=lqr(A,B,params.Q,params.R);
params.K=K0;
%% Data Preparation from Trajectory files
 
 traj_pts=traj_nom(1:end,:);
traj_nom=traj_pts(:,2:11);
 t_nom=traj_pts(:,1);

%% get the trajectory using new gains 
global uu
uu = [params.u0.'];
y0 = params.x0;

t0=0;

nomTraj=[t_nom, traj_nom];
final_time=t_nom(end);

Tf=80; %total time
dt=0.1; 
X=params.x0;
params.wind=wind;
sim_num=1;
data_set=[];
for jj=1:sim_num
    jj
%    wind_idx = (sim_num-1)*10+1;
%    wind_end=wind_idx+Tf/dt;
%    params.wind=wind(wind_idx:wind_end+1,:);
    [t,y]=sdef(Tf,dt,X,params,nomTraj);
    current_data=[t',y];
    data_set=[data_set;current_data];

end

%% input data collect only for single run
u_set=uu(1:4:end,:);
gain_set=kk(1:4:end,:);
data_set=[data_set,u_set];


%% plot trajectory
% 
% plot3(current_data(:,2),current_data(:,3),current_data(:,4))
% hold on
% plot3(traj_nom(:,1),traj_nom(:,2),traj_nom(:,3))
% % 
figure()
plot(current_data(:,1),current_data(:,2))
hold on
% plot(t_nom,traj_nom(:,1))
xlabel('time [s]')
ylabel('x [m]')
title('trajectory in x direction')
legend('Actual','Nominal')
xlim([0,110])


figure()
plot(current_data(:,1),current_data(:,3))
hold on
% plot(t_nom,traj_nom(:,2))
xlabel('time [s]')
ylabel('y [m]')
title('trajectory in y direction')
legend('Actual','Nominal')
xlim([0,110])


figure()
plot(current_data(:,1),current_data(:,4))
hold on
% plot(t_nom,traj_nom(:,3))
xlabel('time [s]')
ylabel('z [m]')
title('trajectory in z direction')
legend('Actual','Nominal')
xlim([0,110])


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
function u = finite_controller(y, params,nomTraj,t,iter)
global uu kk
p = reshape(y(1:3),[3 1]);
q = reshape(y(4:7),[4 1]);
v = reshape(y(8:10),[3 1]);%relative velocity
x=[p;q;v];

if iter==1
    prev_u=uu(end,:);%THIS SHOULD BE CONSTANT FOR ONE SIM so only taking the actual previous
else
    prev_u=uu(end-4,:);
end
% xdes=getref(t,nomTraj);
xdes = params.xref_l;
[A,B] = deriveLinSys_wowind(y,prev_u,params);
K_now = lqr(A,B,params.Q,params.R);
kk=[kk;norm(K_now)];

x(8:10)=x(8:10)+params.meanvel.';%updating to track the ground velocity;

u = params.u0- K_now*(x-xdes');

end

%

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

global uu U

%params
Ix = params.Ix;
Iy = params.Iy;
Iz = params.Iz;
J=[Ix , 0 , 0;
    0 , Iy , 0;
    0 , 0  , Iz];

y;
p = reshape(y(1:3),[3 1]); %position
q = reshape(y(4:7),[4 1]);  %orinetation
v = reshape(y(8:10),[3 1]); % relative velocity
wb = reshape(y(11:13),[3 1]); %omega
q=q/norm(q);


u = finite_controller(y(1:10), params,nomTraj,i,iter);
h_u = u;
% if tc
U = lower_controller(wb,h_u,params);  %low level controller tracks omega

tau = U(1:3);
uu=[uu;u'];

cn = [0;0;U(1)];


% wn = reshape(u(1:3),[3 1]);
% cn = [0;0;u(4)];
% 
% G=  params.G;
% W=  params.W;
% % noise=G*chol(W)*G'*rand(10,1);
% wind=params.wind;
% id=uint32(i*100);
% wind_vel=wind(id+1,:);
% wind_vel(3)=-wind_vel(3);%positive up data

%constant wind
% wind_vel=params.meanvel;
% wind_vel_next=wind_vel;

%LES wind
wind=params.wind;
id=uint32(i*100);
wind_vel_next=wind(id+1,:);

dp=double(v)+wind_vel_next';
dq = double(0.5*quater_mult(q)*[0;wb]);
dvc = double([0;0;params.m*params.g] + quater_rot(q)*cn);
R=quater_rot(q);
dvd=  0.5*params.rho*R*params.D*R'*norm(-v)*-v;
dv=(dvc+dvd)/params.m;
dwb = J\(tau'-cross(wb,(J*wb)));

dx = [dp;dq;dv;dwb];


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

function tau_c = lower_controller(w_c, h_u,P)

  %inputs
  %w_c--- current omega
  %desired-desired(1)-desired wx comes from high level controller
         %-desired(2)-desired wy
         %-desired(3)-desired wz
         %-desired(4)-desired fc
Ix = P.Ix;
Iy = P.Iy;
Iz = P.Iz;
J=[Ix , 0 , 0;
    0 , Iy , 0;
    0 , 0  , Iz];

P_g=diag([100,100,.001]);
des_wb=h_u(1:3);
prev_wb=w_c;

% disp("true")
% disp(true_u(1))
% disp("des")
% disp(des_wb(1))

tau= J*P_g*(des_wb-prev_wb)+ cross(prev_wb,J*prev_wb);

tau_c=[ h_u(4),tau'];

end
function RPM = rpm_controller(f_u,P)
% properties

L = P.l;
b= P.b;
k = P.k;

A = [b b b b;
    0 L*b 0 -L*b;
    L*b 0 -L*b 0;   %from BMET Code 
    k -k k -k];  %

RPM_squared = pinv(A)*f_u'; %pinv-Pseudoinverse

for i=1:4
    if omega_squared(i) < 0
        omega_squared(i) = 0;
    end
end

RPM = sqrt(omega_squared);


end

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

