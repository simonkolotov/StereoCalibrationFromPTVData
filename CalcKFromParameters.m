function K = CalcKFromParameters(varargin)
%K = CalcKFromParameters(f, Cx, Cy)
%
%Calculate the internal calibration matrix K f, Cx and Cy.

if nargin == 3
    f = varargin{1};
    Cx = varargin{2};
    Cy = varargin{3};
    
elseif nargin == 1
    f = varargin{1}(1);
    Cx = varargin{2}(2);
    Cy = varargin{3}(3);
    
else
    disp('Wrong number of input Parameters at CalcKFromParameters');
    keyboard;
end
    

K = [f 0 Cx; 0 f Cy; 0 0 1];
end

