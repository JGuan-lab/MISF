function  [P,S,Q1,F,Q2,Q3_1,Q3_2,err1,err2]=clustering_MISF(X,M,k,k2,alpha,X1,X2,iter,sita)
%%%%The model----------------
% min{P[l],S[l],S} sumTr(P[l]'X[l]L[l]X[l]'P[l])+alpha
% sum||S[l]||^2+sum(\|S[l]-Q1[l]F\|^2+\|M[l]-Q2[l]F\|^2+\|X1-Q3_1F\|^2+\|X2-Q3_2F\|^2)
% You only need to provide the above three inputs.
% Notation:
% X[l] ... (dl x n) gene-cell data 
% n  ... number of cells            
% k ... number of features (the number of cluster)
% k2 ...the number of cluster)
% alpha ... Regularization parameter
% iter ... The maximum number of iterations
%%%

    
    err1 = zeros(iter,1);err2 = zeros(iter,6);
    S1_err = zeros(iter,1);
    S2_err = zeros(iter,1);
    M1_err = zeros(iter,1);
    M2_err =zeros(iter,1);
    X1_err = zeros(iter,1);
    X2_err = zeros(iter,1);

    m = length(X);%2
    [d{1},n] = size(X{1});
    [d{2},n] = size(X{2});
%     P{1}=rand(d{1},k);
%     P{2}=rand(d{2},k);
    [dd{1},nn]=size(X1);
    [dd{2},nn]=size(X2);
    for l=1:m
%         options = [];
%         option.Metric = 'Cosine';%European distance
%         options.NeighborMode = 'KNN';%KNN
%         options.k = 5;%5 nearest neighbors
%         options.WeightMode = 'Cosine';%The weights are 0 or 1, and can be replaced with 'HeatKernel', 'Cosine' 
% 
%         S{l} = constructW(X{l}',options);%%%%%%X is the input matrix
        S{l} = M{l};
        Ds{l} = diag(sum(S{l},2));
	    Ls{l} = Ds{l}-S{l};
        [U,V,D] = svds(S{l},k2);
        Q1{l} = abs(U*sqrt(V));
        fff{l} = abs(sqrt(V)*D');
        [U1,V1,D1] = svds(M{l},k2);
        Q2{l} = abs(U1*sqrt(V1));
        Y = X{l}*Ls{l}*X{l}';
        [Vy,Dy] = eig(Y);
        dy = diag(Dy);
        [t,v] = sort(dy,'ascend');
        P{l} = Vy(:,v(1:k));
    end
    Q3_1=ones(dd{1},k).*0.4;
    Q3_2=ones(dd{2},k).*1;
    F = (fff{1}+fff{2})/2;
for o = 1:iter
%%%%%--------------Update variables P{l} by iteration------------
    for l = 1:m
        %%========P{l}=========
        Y = X{l}*Ls{l}*X{l}';
        [Vy,Dy] = eig(Y);
        dy = diag(Dy);
        [t,v] = sort(dy,'ascend');
        P{l} = Vy(:,v(1:k));
        
    %%%%%====================S{l}======================
    wW = P{l};xX = X{l};
        for i = 1:n
            for j = 1:n
                f(i,j) = norm(wW'*xX(:,i)-wW'*xX(:,j),'fro')^2;
            end
        end
    f=f./2;
    S{l} = S{l}.*((Q1{l}*F)./((alpha+1)*S{l}+f));
    idx=find(isnan(S{l}));
    S{l}(idx)=0;
    S{l} = mapminmax(S{l}, 0, 1);
    end
     %%%%%--------------Update variables B{l},H{l},F by iteration------------
    for l = 1:m
        ssS{l} = (S{l}+S{l}')/2;
       Q1{l} = Q1{l}.*((ssS{l}*F')./(Q1{l}*F*F')); 
       idx3=find(isnan(Q1{l}));
       Q1{l}(idx3)=0;
    end
   
    for l = 1:m
       Q2{l} = Q2{l}.*((M{l}*F')./(Q2{l}*F*F')); 
       idx3=find(isnan(Q2{l}));
       Q2{l}(idx3)=0;
    end
    
    Q3_1 = Q3_1.*((X1*F')./(Q3_1*F*F'));
    Q3_2 = Q3_2.*((X2*F')./(Q3_2*F*F'));
    idx2=find(isnan(Q3_1));
    Q3_1(idx2)=0;
    idx2=find(isnan(Q3_2));
    Q3_2(idx2)=0;

    ff1 = zeros(k2,n); ff2=zeros(k2,n);
    for l = 1:m
        ssS{l} = (S{l}+S{l}')/2;
        ff1 = ff1+Q1{l}'*ssS{l}+Q2{l}'*M{l};
        ff2 = ff2+Q1{l}'*Q1{l}*F+Q2{l}'*Q2{l}*F;
    end
    ff1=ff1+Q3_1'*X1+Q3_2'*X2;
    ff2=ff2+Q3_1'*Q3_1*F+Q3_2'*Q3_2*F;
    F=F.*((ff1)./(ff2));

%%%%%%%%%%%%%%%-------------Error-----------------------
    ee =norm(S{1}-Q1{1}*F,'fro')+norm(S{2}-Q1{2}*F,'fro')+norm(M{1}-Q2{1}*F,'fro')+norm(M{2}-Q2{2}*F,'fro')+norm(X1-Q3_1*F,'fro')+norm(X2-Q3_2*F,'fro');
    err1(o,1)=ee;
    %disp([' Iterations ' num2str(o) ' temp1 ' num2str(ee)]);
    %Calculation of relative error
    ee = norm(S{1}-Q1{1}*F,'fro')/norm(S{1},'fro');
    S1_err(o,1) = ee;
    ee = norm(S{2}-Q1{2}*F,'fro')/norm(S{2},'fro');
    S2_err(o,1) = ee;
    ee = norm(M{1}-Q2{1}*F,'fro')/norm(M{1},'fro');
    M1_err(o,1)= ee;
    ee = norm(M{2}-Q2{2}*F,'fro')/norm(M{2},'fro');
    M2_err(o,1) = ee;
    ee = norm(X1-Q3_1*F,'fro')/norm(X1,'fro');
    X1_err(o,1) = ee;
    ee = norm(X2-Q3_2*F,'fro')/norm(X2,'fro');
    X2_err(o,1) = ee ;
    if o>=2
        e1 = abs(S1_err(o,1)-S1_err(o-1,1));
        e2 = abs(S2_err(o,1)-S2_err(o-1,1));
        e3 = abs(M1_err(o,1)-M1_err(o-1,1));
        e4 = abs(M2_err(o,1)-M2_err(o-1,1));
        e5 = abs(X1_err(o,1)-X1_err(o-1,1));
        e6 =  abs(X2_err(o,1)-X2_err(o-1,1));
    else
        e1 = 1;
        e2 = 1;
        e3 = 1;
        e4 = 1;
        e5 = 1;
        e6 = 1;
    end
    %if err1(o,1) <sita
    if  e1<sita && e2<sita && e3<sita && e4<sita && e5<sita && e6<sita
        err2(:,1) = S1_err;
        err2(:,2) = S2_err;
        err2(:,3) = M1_err;
        err2(:,4) = M2_err;
        err2(:,5) = X1_err;
        err2(:,6) = X2_err;
        break;
    else 
        P = P;
        S = S;
    end
end
  err2(:,1) = S1_err;
  err2(:,2) = S2_err;
  err2(:,3) = M1_err;
  err2(:,4) = M2_err;
  err2(:,5) = X1_err;
  err2(:,6) = X2_err;
 
end

