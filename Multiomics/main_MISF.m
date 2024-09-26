clear all;
close all; 
clc;

%%==============Loading data==============%%
%%%%%%% X1 is scRNA-seq data, rows are genes, columns are cell sample,
%%%%%%% X2 is scATAC-seq data(DNA), rows are genes, columns are cell sample
%%%%%%% Note: scCancer, Pbmc, Kidney data are ATAC; mESC is DNA
%load("scCancer.mat")
load("mESC.mat")
%load("Kidney.mat")
%load("Pbmc.mat")
X1=RNA;
%X2=ATAC;   
X2=DNA;
true_label=real_label;  %Pbmc, Kidney have no real labels, the row should be hidden when applying these two data

%%%%%%% Delete the genes that express all 0 on cells    
for j=1:size(X1,1)
    g=X1(j,:);
    [n1,v1]=find(g~=0);
    shu(j)=length(v1);
end
[v2,n2]=find(shu<=0);
X1(n2,:)=[];
shu=[];
for j=1:size(X2,1)
    g=X2(j,:);
    [n1,v1]=find(g~=0);
    shu(j)=length(v1);
end
[v2,n2]=find(shu<=0);
X2(n2,:)=[];
%%%% Normalize the datas with 0 to 1.
X1=mapminmax(X1,0,1);
X2=mapminmax(X2,0,1);

%%%%%% NMF is utilized to extract features to preprocess the scRNA-seq and scATAC-seq datas
X{1}=X1;X{2}=X2;
%%%%%%r needs to be larger than the number of target clusters
r=15;
for l=1:2
  [W,H]=DR_nmf(X{l},r,110);
  XX{l}=H;
end
XX{1}=mapminmax(XX{1},0,1);
XX{2}=mapminmax(XX{2},0,1);

%%==============Constructing a weight matrix==============%%
%%Preset value before constructing weight matri\
options = [];
option.Metric = 'Cosine';
options.NeighborMode = 'KNN';%KNN
options.k =5;%5 nearest neighbors
options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
W1 = constructW(XX{1}',options);
%%==============Constructing a weight matrix==============%%
%%Preset value before constructing weight matrix
options = [];
option.Metric = 'Cosine';
options.NeighborMode = 'KNN';%KNN
options.k =5;%5 nearest neighbors
options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
W2 = constructW(XX{2}',options);

%%%%%% PMI matrix construction 
A{1}=W1;A{2}=W2;
[M{1}] = PMI(A{1},1);
[M{2}] = PMI(A{2},1);
k=2;k2=2;
alpha=1e-2;
%%%%%%%%%clustering_MISFThe last term of the function input is the threshold sita
[P,S,Q1,F,Q2,Q3_1,Q3_2,err1,err2 ]=clustering_MISF(XX,M,k,k2,alpha,X1,X2,600,1.000000e-3); %%% Call the main function to solve the variables

%%%%%%%%%%%% Clustering cell type label
    for e=1:size(F,2) 
        v=F(:,e);
        ma=max(v);
        [s,t]=find(v==ma);
        prel(e)=s;
    end

%%%%%%%%%%%%%%==================Performance evaluation===============================% %
% External evaluation indicators

result1= ClusteringMeasure_new(true_label,prel); %%%% Accuracy ARI
result2=nmi(true_label,prel);%nmi
result3=Fmeasure(true_label,prel);%F1
result=[result1,result2,result3]
    
