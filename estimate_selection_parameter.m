function [sel mu p]=estimate_selection_parameter(Element_sel,lambda)

sel0=1-10.^(-Element_sel(:,6));
X=[ones(size(Element_sel,1),1),Element_sel(:,1:5)];

sel0_logit=log((sel0+0.0001)./(1.0001-sel0));
sel=sel0_logit;
% sel=zeros(size(sel0_logit,1),1);
k=1;
mu=zeros(6,1);
while k<100
    mu_Old=mu;
    sel_Old=sel;
    mu=pinv(X'*X)*X'*sel;
    sel=1/(1+lambda)*(X*mu+lambda*sel0_logit);
    (norm(sel_Old-sel)/norm(sel))+(norm(mu_Old-mu)/norm(mu));
    k=k+1;
end
sel=exp(sel)./(1+exp(sel));

p=coeff_t_test(sel,X*mu,X,mu);