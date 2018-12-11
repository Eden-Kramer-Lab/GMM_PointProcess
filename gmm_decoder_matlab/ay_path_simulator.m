function [X,Ak,Bk,Qk] = ay_path_simulator(Ak,Bk,Sk,Ns,X0,m_radius)
%% this generates a 2D path in a circle
if nargin==0
    Ak = eye(2,2);
    Sk = 0.2*eye(2,2);
    Bk = [0.06*randn;0.06*randn];
    X0 = [0;0];
    m_radius = 200;
    Ns = 2000;
end
X = zeros(Ns,2);
X(1,:)= mvnrnd(X0,Sk);
for i=2:Ns
    do_loop = 1;
    while do_loop
        Xt = mvnrnd(Ak*X(i-1,:)'+Bk,Sk);
        dr = (Xt-X0')*(Xt-X0')'-m_radius*m_radius;
        if dr<0
            do_loop=0;
            X(i,:)=Xt;      
        end
    end
end
Qk = Sk;
save('Path_X','X','Ak','Bk','Qk');
