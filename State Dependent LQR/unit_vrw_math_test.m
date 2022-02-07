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
w=[wx,wy,wz]';
D= [ dx, 0, 0;
    0, dy, 0;
    0, 0, dz];

new_term=R*cross(R'*v,w);
new_term2=cross(R'*v,w);


dn_dv1=[diff(new_term(1),vx) diff(new_term(1),vy) diff(new_term(1),vz);
    diff(new_term(2),vx) diff(new_term(2),vy) diff(new_term(2),vz);
    diff(new_term(3),vx) diff(new_term(3),vy) diff(new_term(3),vz)];


dn_dv2=[diff(new_term2(1),vx) diff(new_term2(1),vy) diff(new_term2(1),vz);
    diff(new_term2(2),vx) diff(new_term2(2),vy) diff(new_term2(2),vz);
    diff(new_term2(3),vx) diff(new_term2(3),vy) diff(new_term2(3),vz)];%WILL TAKE THIS END THEN MULTIPLY R




%%

dn_dq=[diff(new_term(1),qw)  diff(new_term(1),qx) diff(new_term(1),qy) diff(new_term(1),qz);
       diff(new_term(2),qw)   diff(new_term(2),qx) diff(new_term(2),qy) diff(new_term(2),qz);
       diff(new_term(3),qw)  diff(new_term(3),qx) diff(new_term(3),qy) diff(new_term(3),qz)];
   
   
%%
dn_dw=[diff(new_term2(1),wx) ;
    diff(new_term2(2),wy) ;
    diff(new_term2(3),wz) ];%WILL TAKE THIS END THEN MULTIPLY R