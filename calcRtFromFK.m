function [finalR, finalt, finalF, sigmaD] = calcRtFromFK(F, K, points1, points2)
% [R t sigmaD] = calcRtFromFK(F, K)
%
% Calculate R and t from F given K
% 
% calculation is based on extraction of R and t from E, which is created
% from F by K, and then nonlinear optimization of the results to suit them
% to the points.
%
% F is 3x3 fundamental matrix
%
%       [f 0 Cx
% K is   0 f Cy
%        0 0  1]
%
% points1 and points2 is the list of corresoinding point seen from each camera (nPointsx2), of the form [[x] [y]]
%
%
% R is the rotation part of camera2 relative projection matrix (R' is the relative rotation of cam2 towards cam1)
%
% t is the translation part of camera2 relative projection matrix (-R'*t is the relative translation of cam2)
%
% sigmaD is the sum of the double sided geometric distances over all points
% from their corresponding epipolar lines, as calculated from F.

E = K'*F*K; %E = K2'*F*K1

W = [0 -1 0; 1 0 0; 0 0 1];

[U, S, V] = svd(E);

R_1 = U*W*V'*det(U*W*V');
R_2 = U*W'*V'*det(U*W'*V');

t_1 = U(:,3)*S(1);
t_2 = -t_1;

R = R_1; t = t_1;
H = [R t];
X = triangulate(K*[eye(3) zeros(3,1)], points1, K*H, points2);
x1 = K*[R t]*[X; ones(1,size(X, 2))];
if ~prod(double([X(3,:) x1(3,:)]<0)) %not all z coordinates are < 0
    
    R = R_1; t = t_2;
    H = [R t];
    X = triangulate(K*[eye(3) zeros(3,1)], points1, K*H, points2);
    x1 = K*[R t]*[X; ones(1,size(points1, 2))];
    if ~prod(double([X(3,:) x1(3,:)]<0))
        
        R = R_2; t = t_1;
        H = [R t];
        X = triangulate(K*[eye(3) zeros(3,1)], points1, K*H, points2);
        x1 = K*[R t]*[X; ones(1,size(points1, 2))];
        if ~prod(double([X(3,:) x1(3,:)]<0))
            
            R = R_2; t = t_2;
            H = [R t];
            X = triangulate(K*[eye(3) zeros(3,1)], points1, K*H, points2);
            x1 = K*[R t]*[X; ones(1,size(points1, 2))];
            if ~prod(double([X(3,:) x1(3,:)]<0))
                disp({'Can''t choose [R|t] based on negative z coordinate.'; 'Choosing by rotation and translation constraints'});
                
                if ~prod(double(abs(RotationAngles(R_1'))<pi/2)) %not all 3 rotation angles are under pi/2
                    R = R_2;
                else
                    R = R_1;
                end
                   
                if [1 0 0]*(-R'*t_1)<0 %first coordinate (X axis) of the second camera relative to the first is negative
                    t = t_2;
                else
                    t = t_1;
                end
            end
            
        end
        
    end
end

% nlinSettings = statset('Display', 'Iter');

nlinSettings = statset('Display', 'Off');

initialParameters = [RotationAngles(R) t];
% firstGuessParameters = nlinfit(points, F(:), @CalcFFromParametersForFit, initialParameters, nlinSettings);

finalParameters = nlinfit({K; [points1(1:2,:)' points2(1:2,:)']}, 0, @FitGoal, initialParameters, nlinSettings);

% finalK = CalcKFromParameters(finalParameters(1,1), finalParameters(2,1), finalParameters(3,1));
finalR = RotationMatrix(finalParameters(:,1));
finalt = finalParameters(:,2);

finalF = CalcFFromParameters([KParametersFromK(K) finalParameters]);

sigmaD = median(CalcGeometricDistances(finalF, points1', points2') + CalcGeometricDistances(finalF', points2', points1'));
% % % % 
% % % % RotationAngles(finalR');
% % % % -finalR'*finalt;


end