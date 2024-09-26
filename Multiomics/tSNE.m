clear;
clc;
%%%%处理Pbmc数据的
% 设置字体为 Times New Roman
set(groot, 'defaultAxesFontName','Times New Roman');
set(groot, 'defaultTextFontName','Times New Roman');

% 加载数据
load("scCancer.mat");

%X1 = F;
X1 = RNA;
X2 = ATAC;
labels1 = real_label ;
labels2 =labels1;
%labels2 = readmatrix("Pbmc_MISF_prel.csv");
% t-SNE 参数
numDims = 2; % 设置降维后的维度为2
perplexity = 30; % perplexity 参数，根据需要调整

% t-SNE 降维
Y1 = tsne(X1', 'Algorithm', 'exact', 'NumDimensions', numDims, 'Perplexity', perplexity);
Y2 = tsne(X2', 'Algorithm', 'exact', 'NumDimensions', numDims, 'Perplexity', perplexity);

% 创建颜色映射
uniqueLabels1 = unique(labels1);
numClasses1 = length(uniqueLabels1);
colors1 = colormap(hsv(numClasses1));
% 
uniqueLabels2 = unique(labels2);
numClasses2 = length(uniqueLabels2);
colors2 = colormap(hsv(numClasses2));

% 绘制 t-SNE 降维图
for i = 1:numClasses1
    classIdx1 = labels1 == uniqueLabels1(i);
    scatter(Y1(classIdx1, 1), Y1(classIdx1, 2), 30, colors1(i, :), 'filled');
    hold on;
end

% 添加图例
%legend(arrayfun(@(x) sprintf('Class %d', x), uniqueLabels1, 'UniformOutput', false), 'Location', 'Best');

% 设置标题和坐标轴标签
title('RNA','FontSize',20);
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
hold off;

% 绘制 t-SNE 降维图
figure;
for i = 1:numClasses2
    classIdx2 = labels2 == uniqueLabels2(i);
    scatter(Y2(classIdx2, 1), Y2(classIdx2, 2), 30, colors2(i, :), 'filled');
    hold on;
end

% 添加图例
%legend(arrayfun(@(x) sprintf('Class %d', x), uniqueLabels2, 'UniformOutput', false), 'Location', 'Best');

% 设置标题和坐标轴标签
title('ATAC','FontSize',20);
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
hold off;
