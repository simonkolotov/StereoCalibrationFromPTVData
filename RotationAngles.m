function varargout = RotationAngles(R)
%[thetaX, thetaY, thetaZ] = RotationAngles(R)
%
%or
%
%Angles = RotationAngles(R)
%
%Calculate the three Euler angles from the 3d rotation matrix R, and return
%them in separate variables or as components of a single vector.

thetaX = atan2(R(3,2), R(3,3));

thetaY = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));

thetaZ = atan2(R(2,1), R(1,1));

if (nargout == 3)
    varargout{1} = thetaX;
    varargout{2} = thetaY;
    varargout{3} = thetaZ;

elseif (nargout == 0)||(nargout == 1)
    varargout{1} = [thetaX, thetaY, thetaZ]';

else
    disp ('Wrong number of output parameters');
    keyboard;
end

end