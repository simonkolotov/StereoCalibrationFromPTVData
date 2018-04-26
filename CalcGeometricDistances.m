function distances = CalcGeometricDistances(fundMat, P1, P2)
%geometricDistances = CalcGeometricDistances(fundMat, P1, P2)
%
%Calculates the geometric desitances between epipolar lines
%   calculated from fundMat and P1, and the corresponding points
%   P2. Normalization must be performed outside of this
%   function.

fp = fundMat * P1';
distances = abs(diag( P2*(fp./(repmat(abs(sqrt(fp(1,:).^2+fp(2,:).^2)),3,1)+eps))) );

end