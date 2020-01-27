black_list=[];
[~,i]=ismember('EPAS1',List)
M=size(A,2);
% [~,i]=ismember('EPAS1',List); 
TG=A(i,:);
EG=Element_gene(ismember(Element_gene(:,1),List(i))==1,:);
Element_name_1=unique(EG(cell2mat(EG(:,3))<300000,2));
N=size(Element_name_1,1);

[d f]=ismember(Element_name_1,EG(:,2));
dist_w=exp(-cell2mat(EG(f,3))/300000);

[d f]=ismember(Element_name_1,Element_name);
Ele_act_1=Ele_act(f,:);
sel_1=sel(f);
REid=f;
TF_binding_1=TF_binding(:,f);
gamaPrior_1=gamaPrior(:,find(ismember(gene,List(i))));
TFid=find(abs(gamaPrior_1)>0.01)';
gamaPrior_1=gamaPrior_1(TFid,:);
TF_1=TF(TFid,:);
TFName_1=TFName(TFid,:);
TF_binding_1=TF_binding_1(TFid,:);
isTF=find(ismember(TFName_1,List(i))==1);
 if length(isTF)>0
    TF_1(isTF,:)=[];
    TFName_1(isTF,:)=[];
    TF_binding_1(isTF,:)=[];
    gamaPrior_1(isTF)=[];
     TFid(isTF)=[];
 end
K=size(TF_1,1);
if N*K*sum(sum(Ele_act_1))>0
[d f]=ismember(Element_name_1,EG(:,2));
betaPrior_1=cell2mat(EG(f,4));
%%
try
[alpha,beta,gama,REid,TFid,Z]=alphabetagama2(betaPrior_1,gamaPrior_1,TF_binding_1,TF_1,dist_w,Ele_act_1,TG,N,M,K,delta,sel_1,Element_name_1,TFName_1,TFid,REid); 
catch
    alpha=zeros(size(Ele_act_1,1),1);
    beta=alpha;
    gama=zeros(size(TF_1,1),1);
    act_state=alpha;
    black_list=[black_list;i];
end
act_state=mean(Z,2);
alpha_beta=alpha.*beta;
TFTG{1,i}=[i*ones(length(gama~=0),1) TFid(find(gama~=0))' sign(gama(gama~=0)) gama(gama~=0)];
ElementTG{1,i}=[i*ones(sum(beta~=0),1) REid(find(beta'~=0)) beta(beta~=0) alpha_beta(beta'~=0) act_state(beta'~=0)];
else
TFTG{1,i}=[];
ElementTG{1,i}=[];
end
toc