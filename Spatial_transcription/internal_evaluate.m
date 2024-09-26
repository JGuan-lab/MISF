clc;clear;
%%==============Loading data==============%%
%%%F is the integration matrix of MISF outputs
%%%prel is the prediction result of MISF
load("F.mat")
load("prel.mat")
X=normalize(F);
X=X';
%X =readmatrix('F.csv');
%prel = readmatrix('prel.csv');
%%==============internal evaluate==============%%
chi=evalclusters(X,prel','CalinskiHarabasz');%
CHI=chi.CriterionValues;
dbi=evalclusters(X,prel','DaviesBouldin');%
DBI=dbi.CriterionValues;
sc=evalclusters(X,prel','silhouette');%
SC=sc.CriterionValues;