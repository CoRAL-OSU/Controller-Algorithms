clear all
syms qx qy qz qw c g wx wy wz vx vy vz dx dy dz p N real 

%short cut
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

Qrot1=Q*Q_bar';
R1=Qrot1(2:end,2:end); 

Qrot2=Q'*Q_bar;
R2=Qrot2(2:end,2:end); 

q1=[qw qx qy qz]';
q2=quat_conj(q1);
R_a=quater_rot(q1);
R_a2=quater_rot(q2);

qdot=quater_mult(q1)*[0;wx;wy;wz];
w_quat=[0;wx;wy;wz];
w_quat1=quat_conj(w_quat);

qdot2=quater_mult_bar(w_quat)*q1; %body fixed frame
qdot2_a=quater_mult_bar(q1)*w_quat; %body reference frame
qdot2_b=quater_mult_bar(q1)*w_quat1; %right hand rule

function R = quater_rot(q)

w = q(1);
x = q(2);
y = q(3);
z = q(4);

Q = [w -x -y -z;...
     x w -z y;...
     y z w -x;...
     z -y x w];

Q_bar =[w -x -y -z;...
     x w z -y;...
     y -z w x;...
     z y -x w];
t = Q_bar'*Q; %=Q*Q_bar
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


function Q_bar = quater_mult_bar(q)

w = q(1);
x = q(2);
y = q(3);
z = q(4);

Q_bar =[w -x -y -z;...
     x w z -y;...
     y -z w x;...
     z y -x w];

end

function q_conj=quat_conj(q)
q_conj=[q(1) -q(2) -q(3) -q(4)]';
end

function q_inv=quat_inv(q)
q_conj=[q(1) -q(2) -q(3) -q(4)]';
q_inv=q_conj/(norm(q))^2;
end
