function [] = showPoints(point, varargin)

color = 'r';
MarkerSize = 10;
ShowAxis = 1;

for i = 1 : size(varargin, 2)
    if strcmp(varargin{i}, 'Color')
        color = varargin{i + 1};
    end
    if strcmp(varargin{i}, 'MarkerSize')
        MarkerSize = varargin{i + 1};
    end
    if strcmp(varargin{i}, 'ShowAxis')
        ShowAxis= varargin{i + 1};
    end
end

plot(point(1, :), point(2, :), '.', 'Color', color,  'MarkerSize', MarkerSize)
axis equal

if ShowAxis == 0
    axis off
end
hold off

end