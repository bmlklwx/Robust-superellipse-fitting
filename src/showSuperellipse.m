function [] = showSuperellipse(x, varargin)
color = 'r';
arclength = 0.1;
ShowAxis = 1;

for k = 1 : size(varargin, 2)
    if strcmp(varargin{k}, 'Color')
        color = varargin{k + 1};
    end
    if strcmp(varargin{k}, 'ShowAxis')
        ShowAxis= varargin{k + 1};
    end
    if strcmp(varargin{k}, 'Arclength')
        arclength= varargin{k + 1};
    end
end
x(1) = max(x(1), 0.007);
[point] = uniformSampledSuperellipse(x, arclength, 0);
point(:, end + 1) = point(:, 1);
plot(point(1, :), point(2, :), 'Color', color)
axis equal
if ShowAxis == 0
    axis off
end

end