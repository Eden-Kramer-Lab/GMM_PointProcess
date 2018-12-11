function gmm_m = ay_gmm_mean(post_gmm)

gmm_m = post_gmm(1).w * post_gmm(1).m;
for u=2:length(post_gmm)
    gmm_m = gmm_m + post_gmm(u).w * post_gmm(u).m;
end
