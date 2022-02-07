syms qx qy qz qw c g wx wy wz vx vy vz dx dy dz p N vrx vry vrz vw1 vw2 vw3 real 
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
vr=[vrx,vry,vrz]';
vw=[vw1, vw2, vw3]';
D= [ dx, 0, 0;
    0, dy, 0;
    0, 0, dz];


% 
% 
% drg_term1=norm(vr)*vr;
% 
% drg_dv1=[diff(drg_term1(1),vrx) diff(drg_term1(1),vry) diff(drg_term1(1),vrz)];
% drg_dv2=[diff(drg_term1(2),vrx) diff(drg_term1(2),vry) diff(drg_term1(2),vrz)];
% drg_dv3=[diff(drg_term1(3),vrx) diff(drg_term1(3),vry) diff(drg_term1(3),vrz)];
% 
% drg_dv= R*D*R'*[drg_dv1;
%     drg_dv2;
%     drg_dv3];
% 
% ddrg_dv= R*D*R'*(norm(vr)*eye(3)+(vr*vr')/norm(vr)); %also jacobian of vdot wrt v


%%
% vr=v-vw;
% drg_term2=norm((vr))*(vr);
% 
% 
% drg2_dv1=[diff(drg_term2(1),vx) diff(drg_term2(1),vy) diff(drg_term2(1),vz)];
% drg2_dv2=[diff(drg_term2(2),vx) diff(drg_term2(2),vy) diff(drg_term2(2),vz)];
% drg2_dv3=[diff(drg_term2(3),vx) diff(drg_term2(3),vy) diff(drg_term2(3),vz)];
% 
% drg2_dv=R*D*R'*[drg2_dv1;
%     drg2_dv2;
%     drg2_dv3];
% 
% %% test the equivalency with value
% % RESULT: THEY ALL ARE SAME
%  ddrg_dv=subs(ddrg_dv,[qw,qx,qy,qz,vrx,vry,vrz,dx,dy,dz],[1,0,0,0,-2.70,0,0,0.218,0.218,0.05]);
% drg_dv=subs(drg_dv,[qw,qx,qy,qz,vrx,vry,vrz,dx,dy,dz],[1,0,0,0,-2.70,0,0,0.218,0.218,0.05]);
% drg2_dv=subs(drg2_dv,[qw,qx,qy,qz,vx,vy,vz,vw1,vw2,vw3,dx,dy,dz],[1,0,0,0,0,0,0,2.70,0,0,0.218,0.218,0.05]);


%% lINEARIZING with respect to q
%Linearizing drag component
%wrt q
q=[qw, qx, qy, qz]';
nq = norm(q);
% q*q.'
ddqu_dq = (eye(4) - (nq^-2)*q*q.')*(nq^-1) ;

drg_term=R*D*R'*norm(v)*v; %V has nothing to do with q so ignoring it while diff
ddrg_q1=[diff(drg_term(1),qw) diff(drg_term(1),qx) diff(drg_term(1),qy) diff(drg_term(1),qz)];
ddrg_q2=[diff(drg_term(2),qw) diff(drg_term(2),qx) diff(drg_term(2),qy) diff(drg_term(2),qz)];
ddrg_q3=[diff(drg_term(3),qw) diff(drg_term(3),qx) diff(drg_term(3),qy) diff(drg_term(3),qz)];
ddrg_q= [ddrg_q1;
    ddrg_q2;
    ddrg_q3];%*ddqu_dq; 
