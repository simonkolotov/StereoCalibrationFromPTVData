function R = RotationMatrix(varargin)
%R = RotationMatrix(thetaX, thetaY, thetaZ)
%or
%R = RotationMatrix([thetaX, thetaY, thetaZ])
%
%Calculate the 3d rotation matrix defined by the three Euler angles

if (nargin == 3)
    thetaX = varargin{1};
    thetaY = varargin{2};
    thetaZ = varargin{3};
elseif (nargin == 1)
    thetaX = varargin{1}(1);
    thetaY = varargin{1}(2);
    thetaZ = varargin{1}(3);
else
    disp('wrong number of input parameters');
    keyboard;
end


RX = [1 0 0; 0 cos(thetaX) -sin(thetaX); 0 sin(thetaX) cos(thetaX)];

RY = [cos(thetaY) 0 sin(thetaY); 0 1 0; -sin(thetaY) 0 cos(thetaY)];

RZ = [cos(thetaZ) -sin(thetaZ) 0; sin(thetaZ) cos(thetaZ) 0; 0 0 1];

R = RZ*RY*RX;

end