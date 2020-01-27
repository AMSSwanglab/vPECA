function [sel P_omiga omiga]=estimate_selection(snp,delta,Element_name,Ele_act)
Element_chr_region=regexp(Element_name,'_','split');
Element_chr=cellfun(@(x) x(:,1),Element_chr_region);
M_1 = containers.Map('KeyType','char','ValueType','any');
M_2 = containers.Map('KeyType','char','ValueType','any');
M_3 = containers.Map('KeyType','char','ValueType','any');
M_4 = containers.Map('KeyType','char','ValueType','any');
M_5 = containers.Map('KeyType','char','ValueType','any');
M_cms = containers.Map('KeyType','char','ValueType','any');
temp=cell2mat(snp(:,3:8));
for i=1:length(snp(:,1))
    x=snp(i,1);
    if ~isKey(M_1,x{:})
        M_1(x{:})=temp(i,1);
        M_2(x{:})=temp(i,2);
        M_3(x{:})=temp(i,3);
        M_4(x{:})=temp(i,4);
        M_5(x{:})=temp(i,6);
        M_cms(x{:})=temp(i,5);
    else
        M_1(x{:})=[M_1(x{:}) temp(i,1)];
        M_2(x{:})=[M_2(x{:}) temp(i,2)];
        M_3(x{:})=[M_3(x{:}) temp(i,3)];
        M_4(x{:})=[M_4(x{:}) temp(i,4)];
        M_5(x{:})=[M_5(x{:}) temp(i,6)];
        M_cms(x{:})=[M_cms(x{:}) temp(i,5)];
    end
end
snp_element_max=[cellfun(@max,values(M_1))',cellfun(@max,values(M_2))',cellfun(@max,values(M_3))',...
    cellfun(@max,values(M_4))',cellfun(@max,values(M_5))',cellfun(@max,values(M_cms))'];
snp_element_mean=[cellfun(@mean,values(M_1))',cellfun(@mean,values(M_2))',cellfun(@mean,values(M_3))',...
    cellfun(@mean,values(M_4))',cellfun(@mean,values(M_5))',cellfun(@mean,values(M_cms))'];

% Element_chr_region=regexp(Element_name,'_','split');
Element_matchsnp_name=Element_chr;
for i=1:length(Element_matchsnp_name)
    Element_matchsnp_name{i}=[Element_chr_region{i}{1},':',Element_chr_region{i}{2},'-',Element_chr_region{i}{3}];
end

element_name_dic=keys(M_1);
[a,b]=ismember(element_name_dic,Element_matchsnp_name);
snp_element_max=double(snp_element_max(a,:));
snp_element_mean=double(snp_element_mean(a,:));
element_name_dic=element_name_dic(a);

[a,b]=ismember(Element_matchsnp_name,element_name_dic);
Element_sel_max=zeros(length(Element_name),size(snp_element_max,2));
Element_sel_mean=Element_sel_max;
Element_sel_max(a,:)=snp_element_max(b(a),:);
Element_sel_mean(a,:)=snp_element_mean(b(a),:);
% 
Element_sel=Element_sel_mean;
mu_1=[];
sel_1=[];
for lambda=0:0.01:1
    [sel mu p]=estimate_selection_parameter(Element_sel,lambda);
    mu_1=[mu_1 mu];
    sel_1=[sel_1 sel];
end
[sel mu p]=estimate_selection_parameter(Element_sel,1);
sel(sel<0.95)=0; 

x=ones(5,1);
time_index=[x;2*x;3*x;4*x;5*x];
time_index=[time_index;time_index];
delta=[zeros(25,1);ones(25,1)];

temp=sel*delta';
omiga=zeros(size(Ele_act,1),max(time_index));
eta=omiga;
P_omiga=ones(size(Ele_act,1),max(time_index));
omiga_expand=zeros(size(Ele_act,1),length(time_index));
eta_expand=omiga_expand;
for i=1:size(Ele_act,1)
    if sel(i)~=0
    for j=1:max(time_index)
        a=regstats(Ele_act(i,time_index==j)',[temp(i,time_index==j)'],'linear');
        b=a.beta;
        b(isnan(b))=0;
        omiga(i,j)=b(2);
        eta(i,j)=b(1);
        P_omiga(i,j)=a.tstat.pval(2);
    end
    end
end
P_omiga(isnan(P_omiga))=1;

fid=fopen('./Output//Element_selection.txt','wt');
for i=1:length(Element_name)
fprintf(fid, '%s\t',Element_name{i});
fprintf(fid, '%f\n',sel(i));
end
fid=fopen('./Output/Element_p_omiga.txt','wt');
for i=1:length(Element_name)
fprintf(fid, '%s\t',Element_name{i});
fprintf(fid, '%f\n',P_omiga(i));
end
end