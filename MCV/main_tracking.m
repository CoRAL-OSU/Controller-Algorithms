clc
clear all
close all
%------ Main: Finite Time Minimum Cost Variance Controller ------%
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
global uu tt
tt=0;
params.g = -9.81;
params.Q = diag([10,10,10,1,1,1,1,0.1,0.1,0.1]);
params.R = diag([1,5,5,0.1]);
params.G=G;
params.W=W;
params.gaama=0.7;  %this is highest 0.70
vw_bar=[2.72;1.752;0.0006];
% vw_bar=[0.000001;0;0];
%vw_bar=-[0.001;0.001;0];
params.meanvel=vw_bar;
params.D=[0.218,0,0;
         0, 0.218 , 0;
         0 , 0, 0.004];
     
params.u0 = [0;0;0.05;7]; %initial wz has effect on the trajectory
params.x0 = [0;0;0;1;0;0;0;-params.meanvel]';
params.xref=[0;0;12;1;0;0;0;-params.meanvel]';

params.m=1.035 * 4/6; 
params.rho=1.15;
params.uref=[0;0;0;params.g*params.m];

%% find initial gain value
[A,B]=deriveLinSys(params.xref,params.u0,params);
K0=lqr(A,B,params.Q,params.R);
K=MCV_infinite(A,B,-K0,params);
params.K=K;
%% Data Preparation from Trajectory files


traj_pts=traj_nom(1:end,:);

for i=1:length(traj_pts)
    vw(i,:) = traj_pts(i,9:end)+(-vw_bar');    
end
traj_nom=[traj_pts(:,2:8),  vw];
traj_track=traj_pts(:,2:end);
params.traj_track=traj_track;
%  save("c_traj_nom_lesmran.mat","traj_nom");
%%
t_nom=traj_pts(:,1);
u_nom=control(t_nom,traj_nom,params);

%%
% figure(4)
%  plot(t_nom,traj_nom(:,1))
%  figure(5)
%  plot(t_nom,traj_nom(:,2))
%  
%  figure(6)
%  plot(t_nom,traj_nom(:,3))

%%
nFit = 7;  %Order of polynomial fittig
nState= 10;
nInput= 4;
nTime = length(t_nom);
for i=1:nState
zFit(i,:) = polyfit(t_nom,traj_nom(:,i),nFit);
end
for i=1:nInput
uFit(i,:) = polyfit(t_nom,u_nom(i,:),nFit);
end

%% check fitting

for i=1:length(t_nom)
x_fit(i)= polyval(zFit(1,:),t_nom(i),1);
end
% 
% plot(t_nom,traj_nom(:,1))
% hold on 
% plot(t_nom,x_fit)
%%

%Function handle for getting linearized system dynamics
linSys = @(t)deriveLinTraj(t,zFit,uFit,params);

F =  diag([20,20,20,.1,.1,.1,.1,0.1,0.1,0.1]);  % Terminal cost on state 
tol = 1e-6;  % Accuracy of ricatti propagation
params.F=F;


Soln = trajectoryFMCV(t_nom,linSys,params,tol);

%% see gain propagation
K_val=[];M=[];H=[];E=[];K1=[];K2=[];K3=[];K4=[];K=[];
for i=1:length(t_nom)
    
    gain=Soln(i).K;
    
    K=[K; gain];
    K1=[K1; K(1+(i-1)*nInput,:)];
    K2=[K2; K(2+(i-1)*nInput,:)];
    K3=[K3; K(3+(i-1)*nInput,:)];
    K4=[K4; K(4+(i-1)*nInput,:)];
    gain_val=norm(gain);
    K_val=[K_val;gain_val]; 
    MNow= Soln(i).M;
    M_val=norm(MNow);
    M=[M;M_val];
    HNow= Soln(i).H;
    H_val= norm(HNow);
    H=[H;H_val];
    err=Soln(i).E;
    E=[E;err];
    
end

figure()
plot(t_nom, K_val)
title("Gain Propagation ||K||") 
xlabel("Time,[s]")
ylabel("Gain, ||K||")

figure()
plot(t_nom, M)
title("Riccati Propagation ||M||") 
xlabel("Time,[s]")
ylabel("Riccati, ||M||")

figure()
plot(t_nom, H)
title("Riccati Propagation ||H||") 
xlabel("Time,[s]")
ylabel("Riccati, ||H||")



%% get the trajectory using new gains 
global uu
uu = [params.u0.'];
y0 = params.x0;


t0=0;

nomTraj=traj_pts;%[t_nom, traj_nom];
final_time=t_nom(end);

Tf=110; %total time
dt=0.01; 
X=params.x0;
gain=K;
params.wind=wind;
sim_num=1;
data_set=[];
for jj=1:sim_num
    jj
[t,y]=sdef(Tf,dt,X,params,gain,nomTraj);
current_data=[t',y];
data_set=[data_set;current_data];
end
%%
% %% 

u_set=uu(1:4:end,:);
gain_set = K_val(1:11001,1);
data_set=[data_set,u_set,gain_set];


%% plot trajectory
% 
% plot3(current_data(:,2),current_data(:,3),current_data(:,4))
% hold on
% plot3(traj_nom(:,1),traj_nom(:,2),traj_nom(:,3))
% % 
figure()
plot(current_data(:,1),current_data(:,2))
hold on
plot(t_nom,traj_nom(:,1))
xlabel('time [s]')
ylabel('x [m]')
title('trajectory in x direction')
legend('Actual','Nominal')

figure()
plot(current_data(:,1),current_data(:,3))
hold on
plot(t_nom,traj_nom(:,2))
xlabel('time [s]')
ylabel('y [m]')
title('trajectory in y direction')
legend('Actual','Nominal')

figure()
plot(current_data(:,1),current_data(:,4))
hold on
plot(t_nom,traj_nom(:,3))
xlabel('time [s]')
ylabel('z [m]')
title('trajectory in z direction')
legend('Actual','Nominal')


%%
function u = finite_controller(y, params,gain_i,nomTraj,t)
global uu tt
p = reshape(y(1:3),[3 1]);
q = reshape(y(4:7),[4 1]);
v = reshape(y(8:10),[3 1]);
x=[p;q;v];

Know=gain_i;
xdes=getref(t,nomTraj);

x(8:10)=x(8:10)+params.meanvel;%track the ground velocity
u = params.u0+ Know*(x-xdes');
uu=[uu;u']; tt=[tt;t];
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

function uref=get_u(i,u_nom)

uref=u_nom(i);
end

function dx = quad_dyn(y,params,gain_i,nomTraj,i)

global uu 
y;
p = reshape(y(1:3),[3 1]); %position
q = reshape(y(4:7),[4 1]);  %orinetation
v = reshape(y(8:10),[3 1]); % relative velocity


u = finite_controller(y, params,gain_i,nomTraj,i);

wn = reshape(u(1:3),[3 1]);
cn = [0;0;u(4)];

G=  params.G;
W=  params.W;
% noise=G*chol(W)*G'*rand(10,1);
wind=params.wind;
id=uint32(i*100);
wind_vel=wind(id+1,:);
wind_vel(3)=-wind_vel(3);%positive up data

dp=double(v)+wind_vel';
dq = double(0.5*quater_mult(q)*[0;wn]);
% dv = double([0;0;-params.g] + quater_rot(q)*cn);
dvc = double([0;0;params.m*params.g]+ quater_rot(q)*cn);
R=quater_rot(q);
dvd= -0.5*params.rho*norm(v)*R*params.D*R'*v;
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
function u = control(t,y,params)

uinit=params.u0;
gain=params.K;

u=[];
uprev=uinit;

for i=1:length(t)-1
  
   delx= y(i+1,:)-y(i,:);
   uc= uprev + gain*delx';
   u= [u, uc];
  
end
u=[u, uc]; %for the last point took the previous point

end

function [t,X]=sdef(Tf,dt,X,params,gain,nomTraj)  % Runge Kutta 4
t = 0:dt:Tf;
iter=0;
for i=1:length(t)-1    
    iter=+1;
    gain_i= gain((iter-1)*4+1:iter*4,:);
    X(i,:);
    f1 = quad_dyn(X(i,:),params,gain_i,nomTraj,t(i)) ;
    f2 = quad_dyn(X(i,:)+dt/2*f1,params,gain_i,nomTraj,t(i));
    f3 = quad_dyn(X(i,:)+dt/2*f2,params,gain_i,nomTraj,t(i));
    f4 = quad_dyn(X(i,:)+dt*f3,params,gain_i,nomTraj,t(i));
    X(i+1,:)=X(i,:)+(dt/6)*(f1+2*f2+2*f3+f4);
end
end

