function [alpha beta gama REid TFid Z]=alphabetagama2(betaPrior_1,gamaPrior_1,TF_binding_1,TF_1,dist_w,Ele_act_1,TG,N,M,K,delta,Sel,Element_name_1,TFName_1,TFid,REid)

for ii=1:N
    gBTF(ii,:)=gamaPrior_1'*((TF_binding_1(:,ii)*ones(1,M)).*TF_1);
end
A_beta=diag(dist_w)*[(Ele_act_1).*gBTF];
A_beta=A_beta';

Indices=[repmat([1:5]',5,1);repmat([1:5]',5,1)];
[B,FitInfo]=lasso(A_beta,TG);
mse=zeros(max(Indices),length(FitInfo.Lambda));
for iii=1:max(Indices)
    [B]=lasso(A_beta(Indices~=iii,:),TG(Indices~=iii),'lambda',FitInfo.Lambda);
    yfit = [ones(sum(Indices==iii),1) A_beta(Indices==iii,:)] * [FitInfo.Intercept; B];
    mse(iii,:) = ones(sum(Indices==iii),1)'*(bsxfun(@minus,TG(Indices==iii)',yfit).^2) / sum(Indices==iii);
end
[iii,idxLambdaMinMSE]=min(sum(mse));
x=B(:,idxLambdaMinMSE)~=0;
beta=B(x,idxLambdaMinMSE);

if sum(x)<3
    [B,FitInfo]=lasso(A_beta,TG,'CV',5);
    x=B(:,FitInfo.IndexMinMSE)~=0;
    beta=B(x,FitInfo.IndexMinMSE);
    if sum(x)==0
        [beta,bint,r,rint,stats]=regress(TG',[ones(size(A_beta,1),1) A_beta]);
        beta=beta(2:end);
        if stats(3)<0.05
            x=ones(length(beta),1)>0;
        end
    end
    if sum(x)==1
        x=B(:,FitInfo.IndexMinMSE)~=0;
        beta=regress(TG',[ones(size(A_beta(:,x),1),1) A_beta(:,x)]);
        beta=beta(2:end);
    end
end


REid=REid(x);
Element_name_1=Element_name_1(x);
betaPrior_1=betaPrior_1(x);
TF_binding_1=TF_binding_1(:,x);
dist_w=dist_w(x);
Ele_act_1=Ele_act_1(x,:);
Sel=Sel(x);

clear gBTF
[N,M]=size(Ele_act_1);
K=size(TF_1,1);
for ii=1:N
        BTF{1,ii}=(TF_binding_1(:,ii)*ones(1,M)).*TF_1;
end
A_gama=zeros(K,M);
for i=1:N
    A_gama=A_gama+dist_w(i)*beta(i)*BTF{1,i}.*(ones(K,1)*Ele_act_1(i,:));
end
A_gama=A_gama';

Indices=[repmat([1:5]',5,1);repmat([6:10]',5,1)];
[B,FitInfo]=lasso(A_gama,TG);
mse=zeros(max(Indices),length(FitInfo.Lambda));
for iii=1:max(Indices)
    [B]=lasso(A_gama(Indices~=iii,:),TG(Indices~=iii),'lambda',FitInfo.Lambda);
    yfit = [ones(sum(Indices==iii),1) A_gama(Indices==iii,:)] * [FitInfo.Intercept; B];
    mse(iii,:) = ones(sum(Indices==iii),1)'*(bsxfun(@minus,TG(Indices==iii)',yfit).^2) / sum(Indices==iii);
end
[iii,idxLambdaMinMSE]=min(sum(mse));
x=B(:,idxLambdaMinMSE)~=0;
minMSEModelPredictors = TFName_1(x);
gama=B(x,idxLambdaMinMSE);

TFid=TFid(x);
gamaPrior_1=gamaPrior_1(x);
TF_binding_1=TF_binding_1(x,:);
TF_1=TF_1(x,:);
TFName_1=TFName_1(x,1);

options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
options.MaxIterations = 1000;
x0=double([zeros(length(beta),1);beta;gama]);
fun=@(x) TGfun(x,TG,Ele_act_1,Sel,delta,TF_binding_1,TF_1,dist_w,betaPrior_1,gamaPrior_1);
x=fminunc(fun,x0,options);
if sum(x)~=0
    alpha=x(1:size(Ele_act_1,1));
    beta=x(size(Ele_act_1,1)+1:2*size(Ele_act_1,1));
    gama=x(2*size(Ele_act_1,1)+1:end);
else
    alpha=zeros(size(Ele_act_1,1),1);
end
Z=Ele_act_1+diag(alpha)*Sel*delta';
end


