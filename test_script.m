close all
clear

x_gt = [max(rand * 2, 0.007), ...
       (rand - 0.5) * 5 + 3, (rand - 0.5) * 5 + 3, ...
       rand * 2 * pi, (rand - 0.5) * 5, (rand - 0.5) * 5];
points =  uniformSampledSuperellipse(x_gt, 0.3, 0);
num_point = size(points, 2);
outlier_ratio = 2;
num_out = round(outlier_ratio * num_point);

sigma = mean(eig((points - mean(points, 2)) * (points - mean(points, 2))'/num_point));
outlier = mvnrnd(mean(points, 2)', 2 * sigma * eye(2), num_out)';
points = [points, outlier];

x = EMS2D(points, 'OutlierRatio', 0.8, 'DebugPlot', true);

figure(1)
showPoints(points)
hold on
showSuperellipse(x, 'Color', 'g')
hold off