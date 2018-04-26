function P = triangulate(H1, p1, H2, p2)
% P = triangulate(H1, p1, H2, p2)
%
% TRIANGULATE computes the 3D point location using 2D camera views
%
% H1: camera matrix of the first camera.
% p1: pixel locations matrix of stacked [x, y]' or [x, y, 1]' columns
% H2: camera matrix of the second camera
% p2: pixel locations matrix of stacked [x, y]' or [x, y, 1]' columns
%    (must have the same number of points as points1)
%
% P: the (x, y, z) coordinates of the reconstructed 3D point. Row vector.

nPoints = size(p1,2);

if (size(p2,2) ~= nPoints)
    error('Point vectors are of different Lengths!')
end

P = NaN(4, nPoints);

for iterPoint = 1:nPoints
    %A
    A = [p1(1,iterPoint)*H1(3,:) - H1(1,:); p1(2,iterPoint)*H1(3,:) - H1(2,:); p2(1,iterPoint)*H2(3,:) - H2(1,:); p2(2, iterPoint)*H2(3,:) - H2(2,:)];
    
    [~, ~, V] = svd(A);
    
    P(:,iterPoint) = reshape(V(:,end), 4, 1);
    
end

P = P(1:3,:)./repmat(P(4,:), 3, 1);