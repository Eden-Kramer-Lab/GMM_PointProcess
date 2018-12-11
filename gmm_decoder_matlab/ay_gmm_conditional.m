function [gmm,Ps]=ay_gmm_conditional(Mxy,Cxy,XY,S)
%% Generate gmm
for i=1:size(Mxy,1)
    X(:,i) = mvnpdf(XY,Mxy(i,:),Cxy);
end
Wt  = fitglm(X,S,'linear','Distribution','poisson','Link','identity','Intercept',false);
%Wt  = lassoglm(X,S,'poisson','Link','identity');
Mt  = Wt.Coefficients.Estimate;
ind = find(Mt>1e-2);
%% Find strong ones
while size(X,2)-length(ind)>0
    X   = X(:,ind);
    Mxy = Mxy(ind,:);
    Wt  = fitglm(X,S,'linear','Distribution','poisson','Link','identity','Intercept',false);
    Mt  = Wt.Coefficients.Estimate;
    ind = find(Mt>1e-2);
end
Ps    = glmval(Mt,X,'identity','constant','off');
for i=1:length(Mt)
    gmm(i).W = Mt(i);
    gmm(i).M = Mxy(i,:);
    gmm(i).C = Cxy;
end