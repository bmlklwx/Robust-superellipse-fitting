function [x, p] = EMS2D(point, varargin)

% Written by Weixiao Liu @ JHU, NUS
% Initialized on Feb 24th, 2023, Singapore
% -------------------------------------------------------------------------
% DESCRIPTION: This algorithm solves for the optimal superquadrics (SQ) fitting of a given
%              point cloud. Probabilistic model is adpot to formulate the problem, and
%              thus is roubust enough to tolerate some amount of outliers. The outlier
%              probability is treated as a hidden random variable and is updated via
%              the Bayes' rule (EM algorithm). The parameters of the superquadric is
%              solved iteratively via maximum likelihood estimation.
%
% INPUT: point - point cloud array (3 x N)
%        varargin(optional):
%        'OutlierRatio'        - prior outlier probability [0, 1) (default: 0.1)
%                                we recommend 0 when dealing with clean point cloud.
%        'MaxIterationEM'      - maximum number of EM iterations (default: 20)
%        'ToleranceEM'         - absolute tolerance of EM (default: 1e-3)
%        'RelativeToleranceEM' - relative tolerance of EM (default: 1e-1)
%        'MaxOptiIterations'   - maximum number of optimization iterations per M (default: 2)
%        'Sigma'               - initial sigma^2 (default: 0 - auto generate)
%        'MaxSwitch'           - maximum number of switches allowed (default: 2)
%        'AdaptiveUpperBound'  - introduce adaptive upper bound to restrict the volume of SQ (default: false)
%        'Rescale'             - normalize the input point cloud (default: true)
%
% OUTPUT: x - fitted superquadrics parameters
%         p - outlier probability of the corresponding points
% -------------------------------------------------------------------------
%% Configuration

% parsing varagin
[para] = parseInputArgs(point, varargin{:});

% set EMS parameters
w = para.OutlierRatio;
iterEM_max = para.MaxIterationEM;
iterEM_min = 5;
toleranceEM = para.ToleranceEM;
relative_toleranceEM = para.RelativeToleranceEM;
adaptive_upper = para.AdaptiveUpperBound;

% optimization settings
options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'Display', 'off', 'MaxIterations', para.MaxOptiIterations);

%% Initialization

% translate the coordinate to the center of mass
point = double(point);
t0 = mean(point, 2);
point = point - t0;

% rescale
if para.Rescale == 1
    max_length = max(max(point));
    scale = max_length / 10;
    point = point / scale;
end

% eigen analysis (principal component analysis) for initializing rotation
[EigenVector, ~] = EigenAnalysis(point);
if det(EigenVector) ~= 1
    EigenVector = [EigenVector(:, 2), EigenVector(:, 1)];
end
euler0 = atan2(EigenVector(2,1), EigenVector(1,1));

% initialize scale as median along transformed axis
point_rot0 = EigenVector' * point;
s0 = [median(abs(point_rot0(1, :))), median(abs(point_rot0(2, :)))] * 1;

% initial configuration
x0 = [1, reshape(s0, [1, 2]), euler0, zeros(1, 2)];

%-------------------------------------------
if para.DebugPlot
    figure(1)
    showPoints(point)
    hold on
    showSuperellipse(x0, 'Color', 'g')
    hold off
    title('Init')
    disp('Started in debug mode, please enter in the command window to continue!')
    pause    
end
%-------------------------------------------


% set lower and upper bounds for the superquadrics
upper = 4 * max(max(abs(point)));
lb = [0.0 0.001 0.001 -2*pi -ones(1, 2) * upper];
ub = [2.0 ones(1, 2) * upper  2*pi  ones(1, 2) * upper];

%% EMS Algorithm

% set bounding volume of outlier space
V = (max(point_rot0(1, :)) - min(point_rot0(1, :))) * (max(point_rot0(2, :)) ...
    - min(point_rot0(2, :)));

% outlier probability density
p0 = 1 / V;

% initialize variance for gaussian model
if para.Sigma == 0
    sigma2 = V ^ (1 / 2) / 10;
else
    sigma2 = para.Sigma;
end

% initialize parameters
x = x0;
cost = 0;
switched = 0;
p = ones(1, size(point, 2));

for iterEM = 1 : iterEM_max
    
    % evaluating distance
    [dist] = distance(point, x);
    
    % evaluating corespondence probability
    if w ~= 0
        p = correspendence(dist, sigma2, w, p0);
    end
    
    % adaptive upper bound
    if  adaptive_upper == true
        R_current = angle2rotm(x(4));
        point_rot_current = R_current' * point - R_current' * x(5 : 6)';
        ub_a = 1.1 * [max(abs(point_rot_current(1, :))), ...
            max(abs(point_rot_current(2, :)))];
        ub = [2.0 ub_a  2*pi ub_a];
        lb = [0.001 0.01 0.01 -2*pi -ub_a];
    end
    
    % optimization
    cost_func = @(x) weighted_dist(x, point, p);
    [x_n, cost_n] = lsqnonlin(cost_func, x, lb, ub, options);
    
    % update sigma
    sigma2_n = cost_n / (2 * sum(p));
    
    % evaluate relative cost decrease
    relative_cost = (cost - cost_n) / cost_n;
    if para.DebugPlot
    disp('--------------------------------------------------')
    disp(['relative_cost: ', num2str(relative_cost)])
    disp(['cost: ', num2str(cost)])
    disp(['cost_n: ', num2str(cost_n)])
    end
    %-------------------------------------------
    if para.DebugPlot
        figure(1)
        showPoints(point)
        hold on
        showSuperellipse(x_n, 'Color', 'g')
        hold off
        title(['Trial', num2str(iterEM)])
        pause
    end
    %-------------------------------------------
    
    if (cost_n < toleranceEM && iterEM > 1) || ...
            (relative_cost < relative_toleranceEM ...
            && switched >= para.MaxSwitch && iterEM > iterEM_min)
        x = x_n;
        %-------------------------------------------
        if para.DebugPlot
            figure(1)
            showPoints(point)
            hold on
            showSuperellipse(x, 'Color', 'g')
            hold off
            title('Final')
            pause
        end
        %-------------------------------------------
        break
    end
    
    % set different tolerance for switch and termination
    if relative_cost < relative_toleranceEM && iterEM ~= 1
        if para.DebugPlot
            disp('Enter switching...')
        end
        %-------------------------------------------
        if para.DebugPlot
            figure(1)
            showPoints(point)
            hold on
            showSuperellipse(x_n, 'Color', 'g')
            hold off
            title('Enter switch')
            pause
        end
        %-------------------------------------------
        % activate switching algorithm to avoid local minimum
        switch_success = 0;
        
        % duality similarity
        scale_ratio = x(2) / x(3);
        
        if scale_ratio > 0.6 && scale_ratio < 1.4
            eul_rot = x(4) + 45*pi/180;
            eul_rot = atan(sin(eul_rot)/cos(eul_rot));
            x_candidate = [max(2 - x(1), 1e-2), ((1 - sqrt(2)) * x(1) + sqrt(2)) * min(x(2), x(3)) * ones(1, 2), eul_rot, x(5 : 6)];                        
            %-------------------------------------------
            if para.DebugPlot
                figure(1)
                showPoints(point)
                hold on
                showSuperellipse(x_n, 'Color', 'g')
                showSuperellipse(x_candidate, 'Color', 'b')
                hold off
                title('Candidate in blue')
                pause
            end
            %-------------------------------------------
            % introduce adaptive upper bound
            if  adaptive_upper == 1
                R_current = angle2rotm(x_candidate(4));
                point_rot_current = R_current' * point - R_current' * x_candidate(5 : 6)';
                ub_a = 1.1 * [max(abs(point_rot_current(1, :))), max(abs(point_rot_current(2, :)))];
                ub = [2.0 ub_a  2*pi ub_a];
                lb = [0.0 0.01 0.01 -2*pi -ub_a];
            end
            
            [x_switch, cost_switch] = lsqnonlin(cost_func, x_candidate, lb, ub, options);
            %-------------------------------------------
            if para.DebugPlot
                figure(1)
                showPoints(point)
                hold on
                showSuperellipse(x_n, 'Color', 'g')
                showSuperellipse(x_switch, 'Color', 'b')
                hold off
                title('Candidate updated in blue')
                pause
            end
            %-------------------------------------------
            [dist_switch] = distance(point, x_switch);
            p = correspendence(dist_switch, sigma2, w, p0);
            cost_switch = sum(weighted_dist(x_switch, point, p).^2);
            if cost_switch < min(cost_n, cost)
                x = x_switch;
                cost = cost_switch;
                % update sigma
                sigma2 = cost_switch / (2 * sum(p));
                switch_success = 1;
                %-------------------------------------------
                if para.DebugPlot
                    figure(1)
                    showPoints(point)
                    hold on
                    showSuperellipse(x_n, 'Color', 'g')
                    showSuperellipse(x, 'Color', 'b')
                    hold off
                    title('Candidate switch success')
                    disp('Swith succeed')
                    pause
                end
                %-------------------------------------------
            end
        end
        if switch_success == 0
            cost = cost_n;
            sigma2 = sigma2_n;
            x = x_n;
            %-------------------------------------------
            if para.DebugPlot
                figure(1)
                showPoints(point)
                hold on
                showSuperellipse(x_n, 'Color', 'g')
                hold off
                title('Candidate failed, return to x_n')
                disp('Fail in switch, return to x_n')
                pause
            end
            %-------------------------------------------
        end
        switched = switched + 1;
    else
        cost = cost_n;
        sigma2 = sigma2_n;
        x = x_n;
    end
end
% revert scaling
if para.Rescale == 1
    x(2 : 3) = x(2 : 3) * scale;
    x(5 : 6) = x(5 : 6) * scale;
end

% transform back from the center of mass
x(5 : 6) = x(5 : 6) + t0';
if para.DebugPlot
disp('--------------------------------------------------')
disp('Succeed!')
end
%% Functions
% ------------------eigen analysis-----------------------------------------
    function [EigenVector, EigenValue] = EigenAnalysis(point)
        CovM = point * point' ./ size(point, 2);
        [EigenVector, EigenValue] = eig(CovM);
        EigenVector = flip(EigenVector, 2);
    end

% ------------------distance function--------------------------------------
    function [dist] = distance(X, para)
        % transform pose parameters into R matrix and t vector
        R = [cos(para(4)), -sin(para(4)); sin(para(4)), cos(para(4))];
        t = para(5 : 6);
        % align the point cloud to the superquadrics coordinate
        X_c = R' * X - R' * t';
        % calulate the radial distance of each point
        r_0 = vecnorm(X_c);
        dist = r_0 .* abs((((X_c(1, :) / para(2)) .^ (2)) .^ (1 / para(1)) + ...
            ((X_c(2, :) / para(3)) .^ (2)) .^ (1 / para(1))) .^ (-para(1) / 2) - 1);
    end

% ------------------correspondence calculation ----------------------------
    function [p] = correspendence(dist, sigma2, w, p0)
        c = (2 * pi * sigma2) ^ (- 1);
        const = (w * p0) / (c * (1 - w));
        p = exp(-1 / (2 * sigma2) * dist .^ 2);
        p = p ./ (const + p);
    end

% ------------------weighed distance function -----------------------------
    function [value] = weighted_dist(para, X, p)
        value = p .^ (1 / 2) .* distance(X, para);
        
    end
% ------------------ angle to rotm---------------------
    function [value] = angle2rotm(a)
        value = [cos(a), -sin(a); sin(a), cos(a)];
    end
% ------------------parsing  input arguments-------------------------------
    function [para] = parseInputArgs(point, varargin)
        
        [n_row, n_col] = size(point);
        % check for dimension of input
        if n_row ~= 2
            error('Input point cloud shoud be an array of 2 x N.')
        end
        % check for minumum points
        if n_col < 6
            error('Number of points less than 6.')
        end
        
        % set input parser
        defaults = struct('OutlierRatio', 0.8, ...
            'MaxIterationEM', 50, ...
            'ToleranceEM', 1e-5, ...
            'RelativeToleranceEM', 1e-3, ...
            'MaxOptiIterations', 5, ...
            'Sigma', 0, ...
            'MaxSwitch', 3, ...
            'AdaptiveUpperBound', true, ...
            'Rescale', false,...
            'DebugPlot', false);
        
        parser = inputParser;
        parser.CaseSensitive = false;
        
        parser.addParameter('OutlierRatio', defaults.OutlierRatio, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative' , '>=', 0, '<', 1}));
        parser.addParameter('MaxIterationEM', defaults.MaxIterationEM, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        parser.addParameter('ToleranceEM', defaults.ToleranceEM, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
        parser.addParameter('RelativeToleranceEM', defaults.RelativeToleranceEM, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
        parser.addParameter('MaxOptiIterations', defaults.MaxOptiIterations, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'positive'}));
        parser.addParameter('Sigma', defaults.Sigma, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'nonnegative'}));
        parser.addParameter('MaxSwitch', defaults.MaxSwitch, ...
            @(x)validateattributes(x, {'single', 'double'}, {'real', 'nonsparse', 'nonempty', 'scalar', 'nonnan', 'finite', 'integer', 'nonnegative'}));
        parser.addParameter('AdaptiveUpperBound', defaults.AdaptiveUpperBound, @islogical);
        parser.addParameter('Rescale', defaults.Rescale, @islogical);
        parser.addParameter('DebugPlot', defaults.DebugPlot, @islogical);
        
        parser.parse(varargin{:});
        para = parser.Results;
    end

end