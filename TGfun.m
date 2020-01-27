function [f g]=TGfun(x,TG,Opn_origin,Sel,delta,TF_binding,TF,dist_w,beta0,gama0)
% 
C=1;
mu=1;
% Opn_origin=Ele_act_1;
% Sel=sel_1;
% TF_binding=TF_binding_1;
% TF=TF_1;
% beta0=betaPrior_1;
% gama0=gamaPrior_1;

b=TG'-mean(TG);
alpha=x(1:size(Opn_origin,1));
Opn=Opn_origin+diag(alpha)*Sel*delta';

beta=x(size(Opn,1)+1:2*size(Opn,1));
gama=x(2*size(Opn,1)+1:end);

[N,M]=size(Opn);
K=size(TF,1);
for ii=1:N
        BTF{1,ii}=(TF_binding(:,ii)*ones(1,M)).*TF;
end
for i=1:N
    gBTF(i,:)=gama'*BTF{1,i};
end
A_beta=diag(dist_w)*[Opn.*gBTF];
A_beta=A_beta';
f=1/2*norm(TG'-mean(TG)-A_beta*beta,2)-C*beta0'*beta-C*gama0'*gama+1/2*mu*norm(x);
f=double(f);

A_gama=zeros(K,M);
for i=1:N
    A_gama=A_gama+dist_w(i)*beta(i)*BTF{1,i}.*(ones(K,1)*Opn(i,:));
end
A_gama=A_gama';

b_alpha=b-(diag(dist_w)*(Opn_origin.*gBTF))'*beta;
A_alpha=diag(dist_w.*beta)*((Sel*delta').*gBTF);
A_alpha=A_alpha';

g_alpha=A_alpha'*A_alpha*alpha-A_alpha'*b_alpha+mu*alpha;
g_beta=A_beta'*A_beta*beta-A_beta'*b-C*beta0+mu*beta;
g_gama=A_gama'*A_gama*gama-A_gama'*b-C*gama0+mu*gama;

g=[g_alpha;g_beta;g_gama];
g=double(g);








