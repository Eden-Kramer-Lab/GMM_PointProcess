function [Q,Bs] =ay_gmm_drop_alpha_optimized(P,Q,drop_c)
%% This function drops or merge components of Q by one in reference to Original Mixturte Q
%% Input
%  P: is the original mixture
%  c: is the stopping criteria
%% Output
% Q:  is the merged component
Bs = [];
Q  = P;
%% if there is no c is being defined c~pm * (1-pm)
%% Initialize Algorithm Parameters
% number of P and Q mixtres
n = length(P);
% Likelihood of P over P components mean
p_mp  = zeros(n,n);
pw    = zeros(n,1); 
mps   = zeros(n,length(P(1).m));
for i=1:n
    pw(i)    = P(i).w;
    p_mp(i,:)= pdf_mix(P(i).m,P,n);
    mps(i,:) = P(i).m; 
end
pw     = pw/sum(pw);
qw     = pw;
p_qw   = pw;
p_all  = 0;
q_mq   = p_mp;
p_mq   = p_mp;
q_mp   = p_mp;

%% Main Body, we search for different combination to find the one with lowest pm(1-pm)
drop_info.i = [];
drop_info.Q = 0;
%% two column, first:index second:pf
temp_m     = length(Q);
do_loop = 1;
x_p_mp  = p_mp;
while do_loop
    %% Temporary elemxzSaents
    drop_info.possible = 0;
    ind_m        = 1:temp_m;
    %% Mian loop
    f_max   = -1e6;
    for i =1:temp_m
            %% Build temp_Q,temp_qw
            t_ind  = find(ind_m~=i);
            if (p_all+p_qw(i)) <= (1-drop_c)
                %% Calculate distance (rows are means)
                %% mean of P with Q pdf
                x_q_mp = q_mp(:,t_ind);
                %% mean of Q with P pdf
                x_p_mq = p_mq(t_ind,:);
                %% mean of Q with Q pdf
                x_q_mq = q_mq(t_ind,:);
                x_q_mq = x_q_mq(:,t_ind);
                %% qw
                x_qw   = qw(t_ind)/(1-qw(i));
                %% find the maximum 
                q1 = realmin+(x_q_mp * x_qw);
                q2 = realmin+(x_p_mp * pw);
                t1 = pw'*log(q1./q2);

                q1 = realmin+(x_q_mq * x_qw);
                q2 = realmin+(x_p_mq * pw);
                t2 = - x_qw'*log(q1./q2);
                f  = t1+t2;

                if f >= f_max
                    drop_info.q_mq     = x_q_mq;
                    drop_info.p_mq     = x_p_mq;
                    drop_info.q_mp     = x_q_mp;
                    drop_info.qw       = x_qw;
                    drop_info.t_ind    = t_ind;
                    drop_info.drop_ind = i;
                    f_max   = f;
                    drop_info.possible = 1;
                 end
           end
    end
    
    %% Check whether this is a valid max
    if drop_info.possible
        Bs    = [Bs;f_max];
        p_all = p_all + p_qw(drop_info.drop_ind);
        p_qw  = p_qw(drop_info.t_ind);
        Q     = Q(drop_info.t_ind);
        qw    = drop_info.qw;
        q_mq  = drop_info.q_mq ;
        p_mq  = drop_info.p_mq ;
        q_mp  = drop_info.q_mp;
        temp_m    = length(Q);
        for s=1:temp_m
            Q(s).w=qw(s);
        end
        do_loop = sign(temp_m-1);
    else
        do_loop = 0;
    end
end

function  pf = pdf_mix(x,Mx,dn)
        pf = zeros(dn,1);
        for dd = 1:dn
            pf(dd)= realmin + mvnpdf(x,Mx(dd).m,0.5*(Mx(dd).s+Mx(dd).s'));
        end
end

end
          
          


