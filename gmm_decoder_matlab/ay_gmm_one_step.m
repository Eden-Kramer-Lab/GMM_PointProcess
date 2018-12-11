function post_gmm=ay_gmm_one_step(pre_gmm,Ak,Bk,Sk)
% build one-step gmm
post_gmm = pre_gmm;
for i=1:length(pre_gmm)
    % mean
    post_gmm(i).m = Ak * pre_gmm(i).m + Bk;
    % variance
    post_gmm(i).s = Ak' * pre_gmm(i).s * Ak + Sk;
end
