function tx = SkewMatrix(varargin)
%tx = SkewMatrix(x, y, z)
%
%or
%
%tx = SkewMatrix([x y z])
%
%Calculate the 3x3 skew-symmetric matrix [tx] created from the translation
%vector given by [x y z] or by its coordinates x, y, z.

if nargin == 1
    x = varargin{1}(1);
    y = varargin{1}(2);
    z = varargin{1}(3);

elseif nargin == 3
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};

else
    disp('wrong number of input parameters');
    keyboard;
end

tx = [0 -z y; z 0 -x; -y x 0];

end