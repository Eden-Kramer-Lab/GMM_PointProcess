function [pdf_m,gmm_m,sum_pdf] = ay_gmm_posterior_2d(post_gmm,x_range)


pdf_m = post_gmm(1).w * mvnpdf(x_range,post_gmm(1).m',0.5*(post_gmm(1).s+post_gmm(1).s'));
gmm_m = post_gmm(1).w * post_gmm(1).m;

for u=2:length(post_gmm)
    
    gmm_m = gmm_m + post_gmm(u).w * post_gmm(u).m;
    
    pdf_m = pdf_m + post_gmm(u).w * mvnpdf(x_range,post_gmm(u).m',0.5*(post_gmm(u).s+post_gmm(u).s'));
    
end
sum_pdf = sum(pdf_m);
pdf_m = pdf_m /sum(pdf_m);
