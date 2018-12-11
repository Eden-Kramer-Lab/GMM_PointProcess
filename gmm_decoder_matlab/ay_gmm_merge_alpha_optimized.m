function [Q,Bs] =ay_gmm_merge_alpha_optimized(P,Q,merge_c)
%% This function drops or merge components of Q by one in reference to Original Mixturte Q
%% Input
%  P: is the original mixture
%  c: is the stopping criteria
%% Output
% Q:  is the merged component
Bs = [];
Q  = P;
% cs: is the stopping c 
alpha_range = [linspace(1-2*merge_c,1,10) linspace(max(0,1-10*merge_c),1-2*merge_c-eps,10)];
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
mqs    = mps;
qw     = pw;
q_mq   = p_mp;
p_mq   = p_mp;
q_mp   = p_mp;

%% Main Body, we search for different combination to find the one with lowest pm(1-pm)
merge_info.i = [];
merge_info.j = [];
merge_info.Q = 0;
%% two column, first:index second:pf
temp_m     = length(Q);
do_loop    = 1;
f_max_list = zeros(size(alpha_range));

x_p_mp = p_mp;

while do_loop
    %% Temporary elemxzSaents
    
    merge_info.possible = 0;
    ind_m        = 1:temp_m;
    f_max  = -1e6;
    %% Mian loop
    for i =1:temp_m
        
       for j=i+1:temp_m
            
            %% Build the new mixture
            w1  = qw(i)/(qw(i)+qw(j));
            w2  = qw(j)/(qw(i)+qw(j));
            m0  = w1* Q(i).m + w2* Q(j).m;
            s0  = w1* Q(i).s + w2* Q(j).s + w1*w2*(Q(i).m-Q(j).m)*(Q(i).m-Q(j).m)';
            
            %% Build temp_Q,temp_qw
            t_ind    = [find(ind_m~=i & ind_m~=j) temp_m];
            temp_Q   = Q(t_ind);
            temp_Q(end).m = m0;
            temp_Q(end).s = s0;
            temp_Q(end).w = qw(i)+qw(j);
            
            %% Calculate distance (rows are means)
            %% mean of P with Q pdf
            x_q_mp        = q_mp(:,t_ind);
            x_q_mp(:,end) = mvnpdf(mps,temp_Q(end).m',.5*(temp_Q(end).s+temp_Q(end).s'));
            %% mean of Q with P pdf
            x_p_mq        = p_mq(t_ind,:);
            x_p_mq(end,:) = pdf_mix(temp_Q(end).m,P,n);
            %% mean of Q with Q pdf
            x_q_mq = q_mq(t_ind,:);
            x_q_mq = x_q_mq(:,t_ind);
            x_q_mq(end,:) = pdf_mix(temp_Q(end).m,temp_Q,temp_m-1);
            t_mqs         = mqs(t_ind,:);
            t_mqs(end,:)  = m0;
            x_q_mq(:,end) = mvnpdf(t_mqs,temp_Q(end).m',.5*(temp_Q(end).s+temp_Q(end).s'));
            %% qw
            x_qw          = qw(t_ind);
            x_qw(end)     = qw(i)+qw(j);
            %% find the maximum 
            for s = 1:length(alpha_range)
                 f_max_list(s)=alpha_min(alpha_range(s));
            end
            [~,ind]=max(f_max_list);
            if alpha_range(ind) >= 1-merge_c
                 if f_max_list(ind) > f_max
                    merge_info.q_mq = x_q_mq;
                    merge_info.p_mq = x_p_mq;
                    merge_info.q_mp = x_q_mp;
                    merge_info.qw   = x_qw;
                    merge_info.i = i;
                    merge_info.j = j;
                    merge_info.wab   = qw(i)+qw(j);
                    merge_info.alpha = alpha_range(ind);
                    merge_info.Q     = temp_Q;
                    merge_info.mqs   = t_mqs;
                    f_max  =f_max_list(ind);
                    merge_info.possible = 1;
                    
                 end
            end
        end
    end
    %% Check whether this is a valid max
    if merge_info.possible
        Bs   = [Bs;f_max];
        Q    = merge_info.Q;
        m    = length(Q);
        
        qw   = merge_info.qw;
        wab  = merge_info.wab;
        qw   = qw * (1-merge_info.alpha*wab)/(1-wab);
        qw(end) = merge_info.alpha * wab;
        qw   = qw/sum(qw);
        
        for s=1:m
            Q(s).w=qw(s);
        end
        q_mq = merge_info.q_mq ;
        p_mq = merge_info.p_mq ;
        q_mp = merge_info.q_mp;
        mqs  = merge_info.mqs;
        
        do_loop = sign(m-1);
        temp_m  = m;
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

function  f = alpha_min(x)
          wab  = x_qw(end);
          tqw  = x_qw * (1-x*wab)/(1-wab);
          tqw(end) = x * wab;
          
          q1 = realmin+(x_q_mp * tqw);
          q2 = realmin+(x_p_mp * pw);
          t1 = pw'*log(q1./q2);
        
          q1 = realmin+(x_q_mq * tqw);
          q2 = realmin+(x_p_mq * pw);
          t2 = - tqw'*log(q1./q2);
          f  = t1+t2;
end
    

end
          
          


