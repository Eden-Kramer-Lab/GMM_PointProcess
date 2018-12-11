function post = ay_exact_one_step(pre,Ak,Bk,Sk,Xs)
% build one-step gmm
post = pre;
for i=1:length(Xs)
    % mean
    temp = pdf('normal',Xs(i),Ak*Xs+Bk,sqrt(Sk));
    % variance
    post(i) = pre*temp';
end
post = post/sum(post);
