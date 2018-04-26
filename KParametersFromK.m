function varargout = KParametersFromK(K)
%K = CalcKFromParameters(f, Cx, Cy)
%
%Calculate the internal calibration matrix K f, Cx and Cy.

% K = [f 0 Cx; 0 f Cy; 0 0 1];

f = (K(1,1)+K(2,2))/2; %average of two f

Cx = K(1,3);

Cy = K(2,3);

if (nargout == 3)
    varargout{1} = f;
    varargout{2} = Cx;
    varargout{3} = Cy;

elseif (nargout == 0)||(nargout == 1)
    varargout{1} = [f, Cx, Cy]';

else
    disp ('Wrong number of output parameters');
    keyboard;
end

end

