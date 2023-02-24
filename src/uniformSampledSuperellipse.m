function [point] = uniformSampledSuperellipse(para, arclength ,disp)

sigma = para(1);
scale = para(2 : 3);
xform = para(4 : 6);

threshold = 1e-2;
num_limit = 10000;
theta = zeros(1, num_limit);
seg = zeros(1, num_limit);
theta(1) = 0;
seg(1) = 1;
for m = 2 : num_limit
    [dt, seg_temp] = dtheta(theta(m - 1), arclength, threshold, scale, sigma);
    theta_temp = theta(m - 1) + dt;

    if theta_temp > pi/4
        theta(m) = pi/4;
        break
    else
        if m < num_limit
            theta(m) = theta_temp;
            seg(m) = seg_temp;
        else
            error(['The number of the sampled points exceeds the limit of ', ...
                num2str(num_limit * 4),...
                '. Please increase the arclength or raise the limit'])
        end
    end
end
critical = m + 1;
seg(critical) = 1;
for n = critical + 1 : num_limit
    [dt, seg_temp] = dtheta(theta(n - 1), arclength, threshold, flip(scale), sigma);
    theta_temp = theta(n - 1) + dt;

    if theta_temp > pi/4
        break
    else
        if n < num_limit
            theta(n) = theta_temp;
            seg(n) = seg_temp;
        else
            error(['The number of the sampled points exceeds the limit of ', ...
                num2str(num_limit * 4),...
                '. Please increase the arclength or raise the limit'])
        end
    end
end

num_point = n - 1;
theta = theta(1 : num_point);
seg = seg(1  : num_point);
seg = [seg(1  : critical - 1), flip(seg(critical : end))];

points_fw = angle2points(theta(1 : critical - 1), scale, sigma);
points_bw = flip(angle2points(theta(critical : end), flip(scale), sigma), 2);
point = [points_fw, [points_bw(2, :); points_bw(1, :)]];

point = [point, flip([-point(1, 1 : num_point - 1); point(2, 1 : num_point - 1)], 2), ...
    [-point(1, 2 : end); -point(2, 2 : end)], flip([point(1, 2 : num_point - 1); ...
    -point(2, 2 : num_point - 1)], 2)];

seg = [seg, seg(1 : num_point - 1), seg(2 : end), seg(2 : num_point - 1)];

point = [cos(xform(1)), -sin(xform(1)); sin(xform(1)), cos(xform(1))]...
    * point + [xform(2); xform(3)];

if disp == 1
    figure
    plot(point(1, seg == 1), point(2, seg == 1), '*')
    hold on
    plot(point(1, seg == 2), point(2, seg == 2), '*')
    hold off
    axis equal
end

%------------------ calculation of theta interval----- --------------------
    function [dt, seg] = dtheta(theta, arclength, threshold, scale, sigma)
        if theta < threshold
            dt = abs((arclength / scale(2) + (theta)^(sigma))^(1 / sigma) ...
                - (theta));
            seg = 1;
        else
            dt = arclength / sigma * ((cos(theta) ^ 2 * sin(theta) ^ 2) / ...
                (scale(1) ^ 2 * cos(theta) ^ (2 * sigma) * sin(theta) ^ 4 + ...
                scale(2) ^ 2 * sin(theta) ^ (2 * sigma) * cos(theta) ^ 4))^(1 / 2);
            seg = 2;
        end
    end
%------------------ mapping from angles to points -------------------------
    function [point] = angle2points(theta, scale, sigma)
        point = zeros(2, size(theta, 2));
        point(1, :) = scale(1) .* sign(cos(theta)) .* abs(cos(theta)).^sigma;
        point(2, :) = scale(2) .* sign(sin(theta)) .* abs(sin(theta)).^sigma;
    end

end