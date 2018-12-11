temp = ay_gmm_posterior_1d(pre_gmm,Xs);
plot(Xs,temp/max(temp),'b','linewidth',2);
hold on
temp = ay_gmm_posterior_1d(one_step,Xs);
plot(Xs,temp/max(temp),'k','linewidth',2);

temp = ay_gmm_posterior_1d(px_gmm,Xs);
plot(Xs,temp/max(temp),'r','linewidth',2);

temp = ay_gmm_posterior_1d(post_gmm,Xs);
plot(Xs,temp/max(temp),'y','linewidth',2);

legend('pre','one_step','post','merge')
hold off