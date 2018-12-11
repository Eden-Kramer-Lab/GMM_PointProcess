function Q = ay_gmm_drop(P,Q,drop_c)
%% This function drops or merge components of Q by one in reference to Original Mixturte Q
%% Input
%  P: is the original mixture
%  c: is the stopping criteria
%% Output
% Q:  is the merged component
%% if there is no c is being defined c~pm * (1-pm)
%% Initialize Algorithm Parameters
% number of P and Q mixtres
n = length(P);
% Likelihood of P over P components mean
p_mp  = zeros(n,n);
pw    = zeros(n,1); 
for i=1:n
    pw(i)    = P(i).w;
    p_mp(i,:)= pdf_mix(P(i).m,P,n);
end
ind = find(pw>realmin);
pw  = pw(ind);
Q   = P(ind);
p_mp= p_mp(ind,:);
p_mp= p_mp(:,ind);
pw    = pw./sum(pw);
n     = length(Q);
ind_m = 1:n;
f_max_list = zeros(n,1);
for i=1:n
    t_ind    = find(ind_m~=i);
   
    x_p_mp = p_mp;

    x_q_mq = p_mp(t_ind,:);
    x_q_mq = x_q_mq(:,t_ind);
    
    x_p_mq = p_mp(t_ind,:);
    
    x_q_mp  = p_mp(:,t_ind);
    
    x_qw    = pw(t_ind)/(1-pw(i));
   
    f_max_list(i)=alpha_min();
end
[temp,max_ind] = sort(f_max_list,'descend');
ind_stop   = find(cumsum(pw(max_ind))>= 1-drop_c);
Q          = P(max_ind(ind_stop));



function  pf = pdf_mix(x,Mx,dn)
        pf = zeros(dn,1);
        for dd = 1:dn
            pf(dd)= realmin + mvnpdf(x,Mx(dd).m,0.5*(Mx(dd).s+Mx(dd).s'));
        end
end

function  f = alpha_min()
          
          q1 = realmin + (x_q_mp * x_qw);
          q2 = realmin + (x_p_mp * pw);
          t1 =pw'*log(q1./q2);
        
          q1 = realmin + (x_q_mq * x_qw);
          q2 = realmin + (x_p_mq * pw);
          t2 = - x_qw'*log(q1./q2);
          f  = t1+t2;
end
    

end
          
          


