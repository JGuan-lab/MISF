clear,clc;
load("HPOA.mat");
X=spatial1_1';
% 示例数据，您需要替换为实际的数据
labels1 = type1_1;
labels2 = readmatrix("Cortex_MISF_prel.csv");

% 创建颜色映射
uniqueLabels1 = unique(labels1);
numClasses1 = length(uniqueLabels1);
colors1 = colormap(hsv(numClasses1));

uniqueLabels2 = unique(labels2);
numClasses2 = length(uniqueLabels2);
colors2 = colormap(jet(numClasses2));

% 创建两个子图
%figure;

% 子图1：根据真实标签着色
%subplot(1, 2, 1);
hold on;
for i = 1:numClasses1
    classIdx1 = labels1 == uniqueLabels1(i);
    scatter(X(classIdx1, 1), X(classIdx1, 2), 30, colors1(i, :), 'filled');
end
hold off;

legend(arrayfun(@(x) sprintf('Class %d', x), uniqueLabels1, 'UniformOutput', false), 'Location', 'BestOutside');
legend('boxoff'); % 隐藏图例的边框
title('t-SNE of HPOA Spatial with true label');
xlabel('X');
ylabel('Y');
% 
% % 子图2：根据预测标签着色
% subplot(1, 2, 2);
% hold on;
% for i = 1:numClasses2
%     classIdx2 = labels2 == uniqueLabels2(i);
%     scatter(X(classIdx2, 1), X(classIdx2, 2), 30, colors2(i, :), 'filled');
% end
% hold off;
% 
% % 将图例放在子图外
% legend(arrayfun(@(x) sprintf('Class %d', x), uniqueLabels1, 'UniformOutput', false), 'Location', 'BestOutside');
% legend('boxoff'); % 隐藏图例的边框
% title('t-SNE of Cortex Spatial with Predicted Labels');
% xlabel('X');
% ylabel('Y');
