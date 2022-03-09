
%% TEST RUN FOR PAPER
A = [1.4 , 0.2 , -0.1;
     -0.2 , 0.8, -0.3;
     0.1 , 0.1 , 0.9];
 
B = [0.1 0.8;
     1.1 0.3;
     0.9 0.5;];
 
B1 = [1.2;
    0.1;
    0.2];

Q = eye(3,3);
R = eye(2,2);

params.Q = Q;
params.R = R;
sys.A = A; sys.B=B; sys.B1= B1;
 [K , G]=elqr(sys,params)