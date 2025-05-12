function k = estimate_k_by_svd(X)
% Estimate the optimal number of clusters K using SVD and the Kneedle method
% Input:
%   X - m × n gene expression matrix (rows: cells, columns: genes)
% Output:
%   k - estimated optimal number of clusters

    % Step 1: Z-score normalization
    X = zscore(X);

    % Step 2: SVD decomposition
    [~, S, ~] = svd(X, 'econ');
    singular_values = diag(S);

    % Step 3: Kneedle method to find the elbow point
    n = length(singular_values);
    x = 1:n;
    y = singular_values;

    % Fit a straight line from the first to the last point
    line_start = [x(1), y(1)];
    line_end = [x(end), y(end)];
    line_vec = line_end - line_start;

    % Compute perpendicular distance from each point to the line
    distances = zeros(n, 1);
    for i = 1:n
        point = [x(i), y(i)];
        vec = point - line_start;
        proj = dot(vec, line_vec) / norm(line_vec)^2 * line_vec;
        perp_vec = vec - proj;
        distances(i) = norm(perp_vec);
    end

    % Step 4: Find the point with the maximum distance
    [~, k] = max(distances);

    % Optional: Visualization of singular value spectrum
    figure;
    plot(x, y, 'o-', 'LineWidth', 2);
    hold on;
    plot(k, y(k), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    title(['Estimated Number of Clusters (K = ' num2str(k) ')']);
    xlabel('Component Index');
    ylabel('Singular Value');
    grid on;
    hold off;

end
