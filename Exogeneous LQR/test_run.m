
%% TEST RUN FOR PAPER
% A = [1.4 , 0.2 , -0.1;
%      -0.2 , 0.8, -0.3;
%      0.1 , 0.1 , 0.9];
%  
% B = [0.1 0.8;
%      1.1 0.3;
%      0.9 0.5;];
%  
% B1 = [1.2;
%     0.1;
%     0.2];
clear all

load("model_whole.mat")

A_data=model.A;
A=A_data(1:13,1:13);
B_data=model.B;
B1=B_data(1:13,1:3);
B=B_data(1:13,4:7);

%%

Q = eye(13,13);
R = eye(4,4);

params.Q = Q;
params.R = R;
sys.A = A; sys.B=B; sys.B1= B1;
 [K , G]=elqr(sys,params)
 %%
 
 C=diag([1,1,1,1,1,1,1,1,1,1,1,1,1]);
model=ss(A,B,C,[],Ts=.001);
model_c=d2c(model);
save("model_c.mat","model_c")
A_new=model_c.A; B_new=model_c.B;
k=lqr(A_new,B_new,Q,R);
k1=dlqr(A,B,Q,R);