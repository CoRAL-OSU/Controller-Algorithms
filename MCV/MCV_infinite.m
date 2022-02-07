
function K = MCV_infinite(A,B,K0,P)
%function provide M and H for coupled riccati eqn solution from won et al.
% input A,B,K,Q,R,G, gamma,k 
% k= number of iteration;
% gamma= weight of variance
% A, B state space 
% Q, R gain matrix
% G, W noise 
% K initial control law
% output MH an array that containS M and H
Q = P.Q;
R = P.R;
G = P.G;
W = P.W;
gama = P.gaama;


k=100000;
eps=10^-4;

for i=1:k
    
  
    if i==1
        Kk=K0;
    end
    
    Ak = A+ B*Kk;
    T= Kk'*R*Kk+Q ;

    %Solve Mk
    %Mk=inv(Ak'+Ak)*(-Kk'*R*Kk- Q)
    Mk=lyap(Ak',T);
    %Solve Hk
    %Hk=inv(Ak'+Ak)*4*(-Mk*G*W*G'*Mk)
    Hk=lyap(Ak',4*Mk*G*W*G'*Mk);
    
    Knext=-inv(R)*B'*(Mk+gama*Hk);
    
    
    %del = norm((Knext-Kk)/Kk); %in paper
    del = norm(Knext-Kk);
    
    Kk=Knext;
    if del < eps
      Kopt= Knext; %optimal found
      K=Kopt;
      disp('optimal found with del')
      M=Mk; H=Hk;
      %disp(del)
      break
    
    else
       % disp('optimal not found')
        K=K0; M=Mk; H=Hk;
    end

end

%disp('optimal not found')
end