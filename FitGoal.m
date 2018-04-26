function goal = FitGoal(varargin)
%goal = FitGoal(parameters, points)
%
%
%Function for the nlinfit. It fits the data to the current
%iteration parameters
%
%NOTE: FitGoal doesn't optimize K, but uses it. so the structure to pass on
%to <points> is {K; points}!!!!

parameters = varargin{1};
K = varargin{2}{1};
Points = varargin{2}{2};

% f = 100/.017; Cx = 512; Cy = 640; %width, half width and half height after DS by 4...
% K = CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K

if (nargin == 2)
    K1 = K;
    K2 = K;
else
    K1 = varargin{3};
    K2 = varargin{4};
end

P1 = [Points(:,1:2) ones(size(Points,1),1)];
P2 = [Points(:,3:4) ones(size(Points,1),1)];

F = CalcFFromParameters([KParametersFromK(K) parameters]);

% goal = sum(CalcGeometricDistances(F, P1, P2) + CalcGeometricDistances(F', P2, P1) );
goal = sum(CalcGeometricDistances(F, P1, P2));

end