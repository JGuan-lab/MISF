clear all;
close all; 
clc;
%%==============Loading data==============%%
%%%%%%% X1 is scRNA-seq data, rows are genes, columns are cell sample
%%%%%%% X2 is spatial data, rows are coordinates, columns are cell sample
%%%%%%% Delete the genes that express all 0 on cells
load("Prostate.mat")
%load("HPOA.mat")
%load("Cortex.mat")
%true_label=type1_1;%Note: The Prostate dataset does not have this real label, so you need to omit this line when applying it.
X1=gene1_1;
X2=spatial1_1;
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

X1=mapminmax(X1,0,1);%%%% Normalize the datas with 0 to 1.
X2=mapminmax(X2,0,1);

%%%%%% NMF is utilized to extract features to preprocess the scRNA-seq
X{1}=X1;X{2}=X2;
r=15; %%%r needs to be larger than the number of target clusters
[W,H]=DR_nmf(X{1},r,110);
XX{1}=H;

%initialize W_L (2-by-k), H_L (k-by-n) by kmeans
[j_ids,W_L] = kmeans(X{2}.',r);
H_L = double(bsxfun(@eq, j_ids(:), 1:max(j_ids)));
H_L = H_L.';
W_L = W_L.';
[~,j_ids] = max(H_L);
H_L = confidence_hl(X{2},W_L,H_L);
XX{2}=H_L;

XX{1}=mapminmax(XX{1},0,1);
XX{2}=mapminmax(XX{2},0,1);

%%==============Constructing a weight matrix==============
%%Preset value before constructing weight matrix
options = [];
option.Metric = 'Cosine';
options.NeighborMode = 'KNN';%KNN
options.k =5;%5 nearest neighbors
options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
W1 = constructW(X1',options);
%%==============Constructing a weight matrix==============
%%Preset value before constructing weight matrix
options = [];
option.Metric = 'Cosine';
options.NeighborMode = 'KNN';%KNN
options.k =5;%5 nearest neighbors
options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace withCosine 'HeatKernel', 'Euclidean' 
W2 = constructW(X2',options);

%%%%%% PMI matrix construction 
A{1}=W1;A{2}=W2;
[M{1}] = PMI(A{1},1);
[M{2}] = PMI(A{2},1);
k=7;k2=7;%%%%% k  is the number of feature, k1 is number of cluster%Cortex:k=6; HPOA:k=9; Prostate:k=7
alpha=0.02;
[P,S,Q1,F,Q2,Q3_1,Q3_2,err1,err2]=clustering_MISF(XX,M,k,k2,alpha,X1,X2,1000,1.000000e-6); %%% Call the main function to solve the variables 
%%%%%%%%%%% Clustering cell type label
    for e=1:size(F,2) 
        v=F(:,e);
        ma=max(v);
        [s,t]=find(v==ma);
        prel(e)=s;
    end
%%%%%%%%%%%%%%==================Performance evaluation===============================
% External evaluation indicators
% result1= ClusteringMeasure_new(true_label,prel);%%%% Accuracy ARI
% result2=nmi(true_label,prel);%nmi
% result3=Fmeasure(true_label,prel);%F1
% result=[result1,result2,result3]


%%internal evaluate,
%internal evaluate,"internal_evaluate.m"
X=normalize(F);
X=X';
prel=prel;
chi=evalclusters(X,prel','CalinskiHarabasz');
CHI=chi.CriterionValues
dbi=evalclusters(X,prel','DaviesBouldin');
DBI=dbi.CriterionValues
sc=evalclusters(X,prel','silhouette');
SC=sc.CriterionValues


function [H_L] = confidence_hl(SL,W_L,H_L)
    [~,n] = size(SL);
    for i = 1:n
        dist = pdist2(SL(:,i).', W_L.');
        [dist,idx] = mink(dist,2);
        comp_dist = dist(1) + dist(2);
        H_L(idx(1),i) = dist(2) / comp_dist;
        H_L(idx(2),i) = dist(1) / comp_dist;
    end
end


