%% MCV A and B matrix calculation

syms qx qy qz qw c g wx wy wz vx vy vz dx dy dz p N real 

%Linearizing about [px,py,pz,qw,qx,qy,qz,vx,vy,vz]=[0,0,1,1,0,0,0,0,-2,0]
Q_bar=[qw , -qx, -qy , -qz;
    qx, qw, qz, -qy;
    qy, -qz, qw,qx;
    qz, qy, -qx, qw]; %from paper

Q=[qw , -qx, -qy , -qz;
    qx, qw, -qz, qy;
    qy, qz, qw,-qx;
    qz, -qy, qx, qw];%from paper


Qrot=Q_bar'*Q;%from paper

%this the rotation matrix
R=Qrot(2:end,2:end); %from paper

C=[0 , 0 , c];
G=[0, 0, g]';
v=[vx , vy , vz]';

D= [ dx, 0, 0;
    0, dy, 0;
    0, 0, dz];


% drg_term=R*D*R'*v*norm(v); %in inertial
drg_term=R*D*R'*norm(v)*v; 
gc_term=G- R*C'; 
vdot=gc_term - drg_term;

w=[wx wy wz]';
q=[qw, qx, qy, qz]';
Q_bar_w = [0, -wx, -wy , -wz;
    wx,0,  wz, -wy;
    wy, -wz,0, wx;
    wz, wy, -wx, 0];

qdot=1/2*Q_bar_w*q; % matrix matched with paper refer to eqn 11

simplify(drg_term);
simplify(vdot);

%Linearization
%partical differentation of term with
%respect to each q i.e qw,qx,qy , qz and then multiply with ddqu_dq
%according to paper
nq = norm(q);

ddqu_dq = (eye(4) - (nq^-2)*q*q.')*(nq^-1) ;

ddp_dv = eye([3 3]);
ddq_dq = 0.5*[ 0, -wx, -wy, -wz;...
              wx,   0,  wz, -wy;...
              wy, -wz,   0,  wx;...
              wz,  wy, -wx,   0]*ddqu_dq;
ddq_dw = 0.5*[-qx, -qy, -qz;...
               qw, -qz, -qy;...
               qz,  qw,  qx;...
              -qy,  qx,  qw];



dgc_q1= [diff(gc_term(1),qw) diff(gc_term(1),qx) diff(gc_term(1),qy) diff(gc_term(1),qz)];
dgc_q2=[diff(gc_term(2),qw) diff(gc_term(2),qx) diff(gc_term(2),qy) diff(gc_term(2),qz)];
dgc_q3=[diff(gc_term(3),qw) diff(gc_term(3),qx) diff(gc_term(3),qy) diff(gc_term(3),qz)];
ddvc_dq=[dgc_q1;
    dgc_q2;
    dgc_q3];%*ddqu_dq ;%this matches paper's refer to eqn 22 

%Linearizing drag component
%wrt q
ddrg_q1=[diff(drg_term(1),qw) diff(drg_term(1),qx) diff(drg_term(1),qy) diff(drg_term(1),qz)];
ddrg_q2=[diff(drg_term(2),qw) diff(drg_term(2),qx) diff(drg_term(2),qy) diff(drg_term(2),qz)];
ddrg_q3=[diff(drg_term(3),qw) diff(drg_term(3),qx) diff(drg_term(3),qy) diff(drg_term(3),qz)];
ddrg_q= [ddrg_q1;
    ddrg_q2;
    ddrg_q3]; 
ddrg_dq=ddrg_q;%*ddqu_dq; 
ddrg_dq;
%wrt v
ddvc_dv=zeros(3,1); %v terms are zeros
ddrg_dv= R*D*R'*(norm(v)*eye(3)+(v*v')/norm(v)); %also jacobian of vdot wrt v

%jacobian of vdot wrt q 
ddv_dq=ddvc_dq-ddrg_dq;

ddv_dc = -[            qw*qy + qx*qz;...
                      qy*qz - qw*qx;...
          qw^2 - qx^2 - qy^2 + qz^2];
%  
% ddqu_dq = subs(ddqu_dq,[qw,qx,qy,qz],[1,0 , 0 ,0]);
% ddv_dq= subs(ddv_dq,[qw,qx,qy,qz,vx,vy,vz,dx,dy,dz,c],[1,0,0,0,-0.0001,0,0,0.218,0.218,0.05,10]);
% ddq_dq=subs(ddq_dq,[qw,qx,qy,qz,wx,wy,wz],[1, 0 ,0,0, 0,0,8]);
% 
% %jacobian of vdot wrt v
% ddrg_dv=subs(ddrg_dv,[qw,qx,qy,qz,vx,vy,vz,dx,dy,dz],[1,0,0,0,-0.0001,0,0,0.218,0.218,0.05]);
% ddv_dv = ddvc_dv - ddrg_dv;
% 
% %terms in B 
% ddq_dw=subs(ddq_dw,[qw,qx,qy,qz],[1,0 , 0,0]);
% ddv_dc=subs(ddv_dc,[qw,qx,qy,qz],[1,0 , 0,0]);
% 
% A = [zeros(3),   zeros(3,4), ddp_dv;...
%      zeros(4,3), ddq_dq,     zeros(4,3);...
%      zeros(3) ,  ddv_dq,     ddv_dv];
% B = [zeros(3), zeros(3,1);...
%        ddq_dw, zeros(4,1);...
% 
%        zeros(3),    ddv_dc];
% A=double(A);
% B=double(B);
% save('A_mat0.mat','A')
% save('B_mat0.mat','B')
%% calculation for K with coupled MH at linearized point
% [n,m]=size(A);
% [n1,m1]=size(B);
% Q=diag([10,10,10,1,1,1,1,0.1,0.1,0.1]);
% R=diag([1,5,5,0.1]);
% % G_val=[     0.5434   sqrt(0.1245)   -sqrt(0.0575);
% %     sqrt(0.1245)    0.4197  -sqrt(0.0314);
% %    -sqrt(0.0575)   -sqrt(0.0314)    0.1236;]
% 
% G=zeros(10,3);
% % G(1,1)= 0.2878  ; G(2,2)= 0.1737;G(3,3)=0.0319;
% G(1,1)=  0.5434 ; G(2,2)= 0.4197;G(3,3)=0.1236;
% W=eye(3,3);
% %W=diag([0.5434,0.4197,0.1736]);
% K0=lqr(A,B,Q,R); %initial K0
% gama=1.25;
% [M H K]=LQR_coupled(A,B,Q,R,G,W,-K0,gama);
% -K0
% K
% gainval0=[M;H;K];
% save("zGnw_K_gama1p25.mat","K")
% save("Gnw.mat","G")
% save("G_cov.mat","G")

% save("K_gamap05.mat","K")

%% 