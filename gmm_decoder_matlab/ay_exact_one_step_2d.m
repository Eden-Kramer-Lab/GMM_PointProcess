function post = ay_exact_one_step_2d(pre,Ak,Bk,Sk,Xs)
% build one-step gmm
post = pre;
for i=1:length(Xs)
    % mean
    temp = mvnpdf(Xs(i,:),(Ak*Xs')'+repmat(Bk',size(Xs,1),1),Sk);
    % variance
    post(i) = pre'*temp;
end
post = post/sum(post);
