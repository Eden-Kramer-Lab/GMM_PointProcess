function [pdf_m,gmm_m,sum_pdf] = ay_gmm_posterior_1d(post_gmm,x_range)

pdf_m = post_gmm(1).w * pdf('normal',x_range,post_gmm(1).m,sqrt(post_gmm(1).s));
gmm_m = post_gmm(1).w * post_gmm(1).m;

for u=2:length(post_gmm)
    
    gmm_m = gmm_m + post_gmm(u).w * post_gmm(u).m;
    
    pdf_m = pdf_m + post_gmm(u).w * pdf('normal',x_range,post_gmm(u).m,sqrt(post_gmm(u).s));
    
end
sum_pdf = sum(pdf_m);
pdf_m = pdf_m /sum(pdf_m);
