function varargout = SkewRestore(tx)
%[transX, transY, transZ] = TranslationDistances(tx)
%
%or
%
%translationDistances = TranslationDistances(tx)
%
%extract the translation parameters from the skew-symmetric translation matrix tx, and return
%them in separate variables or as components of a single vector.


transX = tx(6);

transY = tx(7);

transZ = tx(2);

if (nargout == 3)
    varargout{1} = transX;
    
    varargout{2} = transY;
    
    varargout{3} = transZ;
    
elseif (nargout == 0)||(nargout == 1)
    varargout{1} = [transX, transY, transZ]';
    
else
    disp ('Wrong number of output parameters');
    keyboard;
end

end