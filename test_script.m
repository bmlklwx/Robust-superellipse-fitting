close all
clear
addpath('./src')

% random superellipse
x_gt = [max(rand * 2, 0.01), ...
       (rand - 0.5) * 1 + 2, (rand - 0.5) * 1 + 2, ...
       rand * 2 * pi, (rand - 0.5) * 2, (rand - 0.5) * 2];
   
points =  uniformSampledSuperellipse(x_gt, 0.2, 0);
%-------------partial points----------------------
partial_ratio = 0.6; % keep 60% of the points
k = floor(partial_ratio * size(points, 2));
idx = randi(size(points, 2));
distance = vecnorm(points - points(:, idx));
[~, idx_k] = maxk(distance, k);
points = points(:, idx_k);
%-------------------------------------------------

num_point = size(points, 2);
%-------------add outliers 1 = 100%---------------
outlier_ratio = 1;
num_out = round(outlier_ratio * num_point);
sigma = mean(eig((points - mean(points, 2)) * (points - mean(points, 2))'/num_point));
outlier = mvnrnd(mean(points, 2)', 2 * sigma * eye(2), num_out)';
points = [points, outlier];

%--------------add noise--------------------------
noise = mvnrnd([0 0], 0.003 * eye(2), size(points, 2))';
points = points + noise;

tic
x = EMS2D(points, 'OutlierRatio', 0.9, 'DebugPlot', false);
toc

disp('Ground truth is ')
disp(x_gt)
disp('Fitting result is ')
disp(x)

figure(1)
showPoints(points)
hold on
showSuperellipse(x, 'Color', 'g')
hold off