%% prepare data
% HUVEC_FPKM.txt is a gene by sample expression matrix 
% Symbol.txt is the first column of expression matrix from HUVEC_FPKM.txt 
% Element_opn.txt is an element by sample matrix, where each number is a openness score 
% Element_name.txt is the first column of accessibility matrix from Element_opn.txt 
% element_SNP_use.txt is to restore the SNPs locate in elements with their
% selection scores(-log10 p-value). Here the 3rd to 7th columns are Fst, PBS, iHS, XP-EHH,
% and delta DAF
% In Data_prior.mat, the Element_gene file represent the candidate TG and
% RE pairs, their distance and their correlation in public paired samples
% gamaPrior is a TF by TG matrix representing the correlations between them
% in Public data, where TFName and gene are corresponding row and column's
% name
% TF_binding the binding strength derived from motif binding PWM score with
% TFName by Element_name_TFbidnig dimensions

load("./Data/Data_prior.mat")
Symbol=importdata('./Data/Symbol.txt');
A=dlmread('./Data/HUVEC_FPKM.txt','\t',1,1);
List=importdata('./Data/Symbol.txt');
[d f]=ismember(List,Symbol);
d1=ismember(List,gene);
List=List(d.*d1==1);
A=log2(1+A(f(d.*d1==1),:));
[a,b,c]=unique(List);
List=a;
A=A(b,:);
[d f]=ismember(TFName,List);
TF=A(f(d==1),:);

Opn=dlmread('./Data/Element_opn.txt','\t',0,1);
Element_name=importdata('./Data/Element_name.txt');
[a,b,c]=unique(Element_name);
Element_name=Element_name(b);
Opn=Opn(b,:);
d=ismember(Element_gene(:,1),List);
Element_gene=Element_gene(d==1,:);
d=ismember(Element_name,Element_gene(:,2));
Element_name=Element_name(d==1,:);
Element_gene=Element_gene(ismember(Element_gene(:,2),Element_name)==1,:);
Opn=Opn(d==1,:);
Opn1=Opn;
Opn1(Opn1>100)=100;
Ele_act=Opn1;
M=size(A,2);

[a,b]=ismember(Element_name,Element_name_TFbinding);
Element_name_TFbinding=Element_name_TFbinding(b(a),:);
TF_binding=TF_binding(:,b(a));

delta=[zeros(25,1);ones(25,1)];

%% prepare selection data
x=textscan(fopen('./Data/element_SNP_use.txt'),'%s %s %f32 %f32 %f32 %f32 %f32 %f32');
snp=[x{:,1} x{:,2} num2cell(x{:,3}) num2cell(x{:,4}) num2cell(x{:,5}) num2cell(x{:,6}) num2cell(x{:,7}) num2cell(x{:,8})];
[sel P_omiga omiga]=estimate_selection(snp,delta,Element_name,Ele_act);

%% calculate vPECA for each gene

black_list=[];
for i=1:length(List)
tic
gene_name=List(i);
disp(strcat('It is term for gene number:',num2str(i)))
disp(strcat('The gene name is: ',gene_name{:}))
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
end

%%%%%%%%%%
fid=fopen('./Output/TG_TF_regulation.txt','wt');
for i=1:length(List)
if ismember(i,black_list)==0
for k1=1:size(TFTG{1,i})
fprintf(fid, '%s\t',List{TFTG{1,i}(k1,1),1});
fprintf(fid, '%s\t',TFName{TFTG{1,i}(k1,2),1});
fprintf(fid, '%g\t',TFTG{1,i}(k1,3));
fprintf(fid, '%g\n',TFTG{1,i}(k1,4));
end
end
end
fclose(fid);
fid1=fopen('./Output/TG_RE_regulation.txt','wt');
for i=1:length(List)
    if ismember(i,black_list)==0
for k1=1:size(ElementTG{1,i})
fprintf(fid1, '%s\t',List{ElementTG{1,i}(k1,1),1});
fprintf(fid1, '%s\t',Element_name{ElementTG{1,i}(k1,2),1});
fprintf(fid1, '%g\t',ElementTG{1,i}(k1,3));
fprintf(fid1, '%g\t',ElementTG{1,i}(k1,4));
fprintf(fid1, '%g\n',ElementTG{1,i}(k1,5));
end
    end
end
fclose(fid1);
%%%%%
[TFRE(:,1) TFRE(:,2)]=find(TF_binding>0);
ElementTG_corr=corr(A',Opn1');
TFTG_corr=corr(A',TF');
fid2=fopen('./Output/TG_RE_TF_regulation.txt','wt');
for i=1:length(List)
    if ismember(i,black_list)==0
if size(ElementTG{1,i},1)*size(TFTG{1,i},1)>0
d=ismember(TFRE(:,2),ElementTG{1,i}(:,2)).*ismember(TFRE(:,1),TFTG{1,i}(:,2));
TFRE1=TFRE(d==1,:);
for k1=1:size(TFRE1,1)
fprintf(fid2, '%s\t',List{i,1});
fprintf(fid2, '%s\t',Element_name{TFRE1(k1,2),1});
fprintf(fid2, '%s\t',TFName{TFRE1(k1,1),1});
fprintf(fid2, '%g\t',ElementTG_corr(i,TFRE1(k1,2)));
fprintf(fid2, '%g\n',TFTG_corr(i,TFRE1(k1,1)));
end
end
end
end
fclose(fid2);
