function [A,B] = deriveLinSys(y,u,params)
%GETLINSYS
%    [A,B] = GETLINSYS(X,Y,U)


p = reshape(y(1:3),[3 1]);
q = reshape(y(4:7),[4 1]);
v = reshape(y(8:10),[3 1]);
nq = norm(q);
% q*q.'
ddqu_dq = (eye(4) - (nq^-2)*q*q.')*(nq^-1) ;

wx = u(1);
wy = u(2);
wz = u(3);
c = u(4);

qw = q(1);
qx = q(2);
qy = q(3);
qz = q(4);

ddp_dv = eye([3 3]);
ddq_dq = 0.5*[ 0, -wx, -wy, -wz;...
              wx,   0,  wz, -wy;...
              wy, -wz,   0,  wx;...
              wz,  wy, -wx,   0]*ddqu_dq;
ddq_dw = 0.5*[-qx, -qy, -qz;...
               qw, -qz, -qy;...
               qz,  qw,  qx;...
              -qy,  qx,  qw];
ddc_dq = 2*c*[ qy,  qz,  qw, qx;...
              -qx, -qw,  qz, qy;...
               qw, -qx, -qy, qz]*ddqu_dq;
ddv_dc = [            qw*qy + qx*qz;...
                      qy*qz - qw*qx;...
          qw^2 - qx^2 - qy^2 + qz^2]/params.m;
      
R=quater_rot(q);
D=params.D;
ddrg_v=( -0.5*params.rho*R*D*R'*(norm(v)*eye(3)+(v*v')/norm(v)))/params.m; %jacobian of vdot wrt v     
vx=v(1);
vy=v(2);
vz=v(3);
dx=D(1);
dy=D(2);
dz=D(3);

%updated drag diff with respect to q doesnot inlude normalization term ddqu_dq
ddrg_q=0.5*params.rho*[vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dy*qz*(2*qw*qz - 2*qx*qy) + 4*dz*qy*(2*qw*qy + 2*qx*qz) + 4*dx*qw*(qw^2 + qx^2 - qy^2 - qz^2)) - vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qw*(2*qw*qy - 2*qx*qz) + 2*dy*qx*(2*qw*qz - 2*qx*qy) - 2*dz*qw*(2*qw*qy + 2*qx*qz) + 2*dy*qz*(2*qw*qx + 2*qy*qz) + 2*dx*qy*(qw^2 + qx^2 - qy^2 - qz^2) - 2*dz*qy*(qw^2 - qx^2 - qy^2 + qz^2)) - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dy*qw*(2*qw*qz - 2*qx*qy) - 2*dx*qw*(2*qw*qz + 2*qx*qy) + 2*dz*qx*(2*qw*qy + 2*qx*qz) + 2*dz*qy*(2*qw*qx - 2*qy*qz) - 2*dx*qz*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dy*qz*(qw^2 - qx^2 + qy^2 - qz^2)),   vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qx*(2*qw*qz + 2*qx*qy) + 2*dy*qx*(2*qw*qz - 2*qx*qy) - 2*dz*qw*(2*qw*qy + 2*qx*qz) - 2*dz*qz*(2*qw*qx - 2*qy*qz) + 2*dx*qy*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dy*qy*(qw^2 - qx^2 + qy^2 - qz^2)) - vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qx*(2*qw*qy - 2*qx*qz) + 2*dy*qw*(2*qw*qz - 2*qx*qy) - 2*dy*qy*(2*qw*qx + 2*qy*qz) + 2*dz*qx*(2*qw*qy + 2*qx*qz) - 2*dx*qz*(qw^2 + qx^2 - qy^2 - qz^2) - 2*dz*qz*(qw^2 - qx^2 - qy^2 + qz^2)) + vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dz*qz*(2*qw*qy + 2*qx*qz) - 4*dy*qy*(2*qw*qz - 2*qx*qy) + 4*dx*qx*(qw^2 + qx^2 - qy^2 - qz^2)), vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qy - 2*qx*qz) + 2*dy*qx*(2*qw*qx + 2*qy*qz) - 2*dy*qz*(2*qw*qz - 2*qx*qy) - 2*dz*qy*(2*qw*qy + 2*qx*qz) - 2*dx*qw*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dz*qw*(qw^2 - qx^2 - qy^2 + qz^2)) - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qz + 2*qx*qy) + 2*dz*qw*(2*qw*qx - 2*qy*qz) + 2*dy*qy*(2*qw*qz - 2*qx*qy) - 2*dz*qz*(2*qw*qy + 2*qx*qz) - 2*dx*qx*(qw^2 + qx^2 - qy^2 - qz^2) - 2*dy*qx*(qw^2 - qx^2 + qy^2 - qz^2)) - vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dy*qx*(2*qw*qz - 2*qx*qy) - 4*dz*qw*(2*qw*qy + 2*qx*qz) + 4*dx*qy*(qw^2 + qx^2 - qy^2 - qz^2)),   vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qz*(2*qw*qy - 2*qx*qz) - 2*dy*qw*(2*qw*qx + 2*qy*qz) - 2*dy*qy*(2*qw*qz - 2*qx*qy) + 2*dz*qz*(2*qw*qy + 2*qx*qz) + 2*dx*qx*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dz*qx*(qw^2 - qx^2 - qy^2 + qz^2)) - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qz*(2*qw*qz + 2*qx*qy) + 2*dz*qx*(2*qw*qx - 2*qy*qz) - 2*dy*qz*(2*qw*qz - 2*qx*qy) - 2*dz*qy*(2*qw*qy + 2*qx*qz) - 2*dx*qw*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dy*qw*(qw^2 - qx^2 + qy^2 - qz^2)) + vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dy*qw*(2*qw*qz - 2*qx*qy) + 4*dz*qx*(2*qw*qy + 2*qx*qz) - 4*dx*qz*(qw^2 + qx^2 - qy^2 - qz^2));
vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dx*qz*(2*qw*qz + 2*qx*qy) + 4*dz*qx*(2*qw*qx - 2*qy*qz) + 4*dy*qw*(qw^2 - qx^2 + qy^2 - qz^2)) - vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qz + 2*qx*qy) - 2*dy*qw*(2*qw*qx + 2*qy*qz) + 2*dz*qw*(2*qw*qx - 2*qy*qz) + 2*dx*qz*(2*qw*qy - 2*qx*qz) - 2*dy*qx*(qw^2 - qx^2 + qy^2 - qz^2) + 2*dz*qx*(qw^2 - qx^2 - qy^2 + qz^2)) - vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dy*qw*(2*qw*qz - 2*qx*qy) - 2*dx*qw*(2*qw*qz + 2*qx*qy) + 2*dz*qx*(2*qw*qy + 2*qx*qz) + 2*dz*qy*(2*qw*qx - 2*qy*qz) - 2*dx*qz*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dy*qz*(qw^2 - qx^2 + qy^2 - qz^2)),   vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qx*(2*qw*qz + 2*qx*qy) + 2*dy*qx*(2*qw*qz - 2*qx*qy) - 2*dz*qw*(2*qw*qy + 2*qx*qz) - 2*dz*qz*(2*qw*qx - 2*qy*qz) + 2*dx*qy*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dy*qy*(qw^2 - qx^2 + qy^2 - qz^2)) - vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qy - 2*qx*qz) + 2*dy*qx*(2*qw*qx + 2*qy*qz) - 2*dx*qz*(2*qw*qz + 2*qx*qy) - 2*dz*qx*(2*qw*qx - 2*qy*qz) - 2*dy*qw*(qw^2 - qx^2 + qy^2 - qz^2) + 2*dz*qw*(qw^2 - qx^2 - qy^2 + qz^2)) + vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dx*qy*(2*qw*qz + 2*qx*qy) + 4*dz*qw*(2*qw*qx - 2*qy*qz) - 4*dy*qx*(qw^2 - qx^2 + qy^2 - qz^2)), vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dy*qy*(2*qw*qx + 2*qy*qz) - 2*dx*qx*(2*qw*qy - 2*qx*qz) - 2*dx*qw*(2*qw*qz + 2*qx*qy) + 2*dz*qy*(2*qw*qx - 2*qy*qz) + 2*dy*qz*(qw^2 - qx^2 + qy^2 - qz^2) + 2*dz*qz*(qw^2 - qx^2 - qy^2 + qz^2)) - vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qz + 2*qx*qy) + 2*dz*qw*(2*qw*qx - 2*qy*qz) + 2*dy*qy*(2*qw*qz - 2*qx*qy) - 2*dz*qz*(2*qw*qy + 2*qx*qz) - 2*dx*qx*(qw^2 + qx^2 - qy^2 - qz^2) - 2*dy*qx*(qw^2 - qx^2 + qy^2 - qz^2)) + vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dx*qx*(2*qw*qz + 2*qx*qy) - 4*dz*qz*(2*qw*qx - 2*qy*qz) + 4*dy*qy*(qw^2 - qx^2 + qy^2 - qz^2)), - vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qz*(2*qw*qz + 2*qx*qy) + 2*dz*qx*(2*qw*qx - 2*qy*qz) - 2*dy*qz*(2*qw*qz - 2*qx*qy) - 2*dz*qy*(2*qw*qy + 2*qx*qz) - 2*dx*qw*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dy*qw*(qw^2 - qx^2 + qy^2 - qz^2)) - vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qw*(2*qw*qy - 2*qx*qz) - 2*dx*qx*(2*qw*qz + 2*qx*qy) + 2*dy*qz*(2*qw*qx + 2*qy*qz) + 2*dz*qz*(2*qw*qx - 2*qy*qz) - 2*dy*qy*(qw^2 - qx^2 + qy^2 - qz^2) - 2*dz*qy*(qw^2 - qx^2 - qy^2 + qz^2)) - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dz*qy*(2*qw*qx - 2*qy*qz) - 4*dx*qw*(2*qw*qz + 2*qx*qy) + 4*dy*qz*(qw^2 - qx^2 + qy^2 - qz^2));
vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dx*qy*(2*qw*qy - 2*qx*qz) + 4*dy*qx*(2*qw*qx + 2*qy*qz) + 4*dz*qw*(qw^2 - qx^2 - qy^2 + qz^2)) - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qz + 2*qx*qy) - 2*dy*qw*(2*qw*qx + 2*qy*qz) + 2*dz*qw*(2*qw*qx - 2*qy*qz) + 2*dx*qz*(2*qw*qy - 2*qx*qz) - 2*dy*qx*(qw^2 - qx^2 + qy^2 - qz^2) + 2*dz*qx*(qw^2 - qx^2 - qy^2 + qz^2)) - vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qw*(2*qw*qy - 2*qx*qz) + 2*dy*qx*(2*qw*qz - 2*qx*qy) - 2*dz*qw*(2*qw*qy + 2*qx*qz) + 2*dy*qz*(2*qw*qx + 2*qy*qz) + 2*dx*qy*(qw^2 + qx^2 - qy^2 - qz^2) - 2*dz*qy*(qw^2 - qx^2 - qy^2 + qz^2)), - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qy - 2*qx*qz) + 2*dy*qx*(2*qw*qx + 2*qy*qz) - 2*dx*qz*(2*qw*qz + 2*qx*qy) - 2*dz*qx*(2*qw*qx - 2*qy*qz) - 2*dy*qw*(qw^2 - qx^2 + qy^2 - qz^2) + 2*dz*qw*(qw^2 - qx^2 - qy^2 + qz^2)) - vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qx*(2*qw*qy - 2*qx*qz) + 2*dy*qw*(2*qw*qz - 2*qx*qy) - 2*dy*qy*(2*qw*qx + 2*qy*qz) + 2*dz*qx*(2*qw*qy + 2*qx*qz) - 2*dx*qz*(qw^2 + qx^2 - qy^2 - qz^2) - 2*dz*qz*(qw^2 - qx^2 - qy^2 + qz^2)) - vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dx*qz*(2*qw*qy - 2*qx*qz) - 4*dy*qw*(2*qw*qx + 2*qy*qz) + 4*dz*qx*(qw^2 - qx^2 - qy^2 + qz^2)), vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qy*(2*qw*qy - 2*qx*qz) + 2*dy*qx*(2*qw*qx + 2*qy*qz) - 2*dy*qz*(2*qw*qz - 2*qx*qy) - 2*dz*qy*(2*qw*qy + 2*qx*qz) - 2*dx*qw*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dz*qw*(qw^2 - qx^2 - qy^2 + qz^2)) + vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dy*qy*(2*qw*qx + 2*qy*qz) - 2*dx*qx*(2*qw*qy - 2*qx*qz) - 2*dx*qw*(2*qw*qz + 2*qx*qy) + 2*dz*qy*(2*qw*qx - 2*qy*qz) + 2*dy*qz*(qw^2 - qx^2 + qy^2 - qz^2) + 2*dz*qz*(qw^2 - qx^2 - qy^2 + qz^2)) + vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dx*qw*(2*qw*qy - 2*qx*qz) + 4*dy*qz*(2*qw*qx + 2*qy*qz) - 4*dz*qy*(qw^2 - qx^2 - qy^2 + qz^2)),   vx*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qz*(2*qw*qy - 2*qx*qz) - 2*dy*qw*(2*qw*qx + 2*qy*qz) - 2*dy*qy*(2*qw*qz - 2*qx*qy) + 2*dz*qz*(2*qw*qy + 2*qx*qz) + 2*dx*qx*(qw^2 + qx^2 - qy^2 - qz^2) + 2*dz*qx*(qw^2 - qx^2 - qy^2 + qz^2)) - vy*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(2*dx*qw*(2*qw*qy - 2*qx*qz) - 2*dx*qx*(2*qw*qz + 2*qx*qy) + 2*dy*qz*(2*qw*qx + 2*qy*qz) + 2*dz*qz*(2*qw*qx - 2*qy*qz) - 2*dy*qy*(qw^2 - qx^2 + qy^2 - qz^2) - 2*dz*qy*(qw^2 - qx^2 - qy^2 + qz^2)) + vz*(abs(vx)^2 + abs(vy)^2 + abs(vz)^2)^(1/2)*(4*dy*qy*(2*qw*qx + 2*qy*qz) - 4*dx*qx*(2*qw*qy - 2*qx*qz) + 4*dz*qz*(qw^2 - qx^2 - qy^2 + qz^2))];
 
ddv_dq=(ddc_dq-ddrg_q*ddqu_dq)/params.m;

%ensure only v component get divided by mass

A = [zeros(3),   zeros(3,4), ddp_dv;...
     zeros(4,3), ddq_dq,     zeros(4,3);...
     zeros(3) ,  ddv_dq,     ddrg_v];
B = [zeros(3), zeros(3,1);...
       ddq_dw, zeros(4,1);...
       zeros(3),    ddv_dc];
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
