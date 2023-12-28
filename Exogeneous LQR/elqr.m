%for exogeneous input test code from the paper

function [K , G]=elqr(sys,params)

A = sys.A;
B = sys.B;
B1 = sys.B1;
Q = params.Q;
R = params.R;

P = idare(A,B,Q,R,[],[])
K = inv(R+B'*P*B)*B'*P*A;
S = inv(A*inv(P-Q)-inv(P))*B1;
G = K*inv(P-Q)*S;
end