classdef CorrespondingPointList < handle
    %CorrespondingPointList is a class defining the List holding the Corresponding Points extracted from SuperTrackList, for the PTV Point Matching Algorithm
    %
    %Data Structure: An array of 1xCamVewGroups (see below). each cell of the array contains a (2xnCamsPerGroup) x nMatchingPoints array of point coordinates,
    %   with [x y, x y, x y...] structure.
    %   **Current Default Group is a Pair**
    %
    %Parameters:
    %     normalizeCoordinates (Default 1) - whether to perform coordinate
    %       normalization when calculating the fundamental matrix
    %
    %     calcFMethod (Default 'IRLS') - the  method of calculation of the
    %       fundamental matrix (see FUNDAMENTAL MATRIX CALCULATION below)
    %
    %     pointWeights - The weights of the points, as returned by the IRLS
    %       method of calculation for the fundamental matrix, and used
    %       later on by the non linear fit. The weight of each point
    %       reflects the confidence of the point pair being corresponding,
    %       i.e. complying to the epipolar constraint
    %
    %CONSTRUCTORS:
    %   CorrespondingPointListInstance = CorrespondingPointList() - Create empty object of CorrespondingPointList Class
    %
    %   CorrespondingPointListInstance = CorrespondingPointList(CPLInstanceOld) - Initialize from another CPL Instance <CPLInstanceOld>
    %
    %   CorrespondingPointListInstance = CorrespondingPointList(superTrackListIn, superSpaceTimeVolumeIn, frameGroup) -
    %       Initialize by extracting points from STL <superTrackListIn> (assuming STL is already matched),
    %       using data from SSTV <superSpaceTimeVolumeIn> in frames defined by <frameGroup>, which is a vector of frame indices (e.g. [1, 15]).
    %
    %METHODS:
    %   REMOVE DOUBLE POINTS
    %       CPL.RemoveDoublePoints() - Remove points occuring more than once in the list
    %
    %   ADD POINTS FROM ANOTHER SLICE
    %       CPL.AddPointsFromAnotherSlice(superTrackListIn, superSpaceTimeVolumeIn, frameGroup)
    %           Add points extracted from STL <superTrackListIn> (assuming STL is already matched),
    %           using data from SSTV <superSpaceTimeVolumeIn> in frames defined by <frameGroup>, which is a vector of frame indices (e.g. [1, 15]).
    %
    %   FUNDAMENTAL MATRIX CALCULATION
    %
    %       fundMat = CPL.CalcFundMat() - Calculate the Fundamental Matrix from the List of Matching Points using the method set in <calcFMethod>.
    %
    %       fundMat = CPL.CalcFundMat(method) - Calculate the Fundamental Matrix from the List of Matching Points, using the method <method>,
    %           which can be 'LS' for simple Least Squares, 'RANSAC', or 'IRLS' for Iteratively Reweighted Least Squares (default).
    %
    %       fundMat = CPL.CalcFundMat(method, normalize) - Same as above, <normalize> is a binary variable stating whether to perform coordinate normalization
    %           or not. (default = 1)
    %
    %       In case of more than two cameras (CPL is an array, and the 4th dimention of the SSTV is of the appropriate size - for pairs, one camera more than the number of
    %           groups), fundMat is a 3D array, having a 3x3 layer for each camera group.
    %
    % % % %     %   point = STL.GetPoint(PointNum, cameraNum, direction) - Retrieve point ([x y]) number <pointNum> in <cameraNum> camera view.
    % % % %     %       In case of more than two Camera Views, <direction> is used to describe the camera position in the pair. <direction> can be 'forward' (default) or
    % % % %     %       'backward'.
    %
    %Simon Kolotov, Ver 2.0, Spring 2013
    
    properties (SetAccess = public, GetAccess = public)
        data
        
        normalizeCoordinates = 1;
        
        calcFMethod = 'IRLS';
        
        pointWeights
        
    end %public Properties
    
    methods
        %CONSTRUCTORS
        function CorrespondingPointListInstance = CorrespondingPointList(varargin)
            %   CorrespondingPointListInstance = CorrespondingPointList() - Create empty object of CorrespondingPointList Class
            %
            %   CorrespondingPointListInstance = CorrespondingPointList(CPLInstanceOld) - Initialize from another CPL Instance <CPLInstanceOld>
            %
            %   CorrespondingPointListInstance = CorrespondingPointList(superTrackListIn, superSpaceTimeVolumeIn, frameGroup) -
            %       Initialize by extracting points from STL <superTrackListIn> (assuming STL is already matched),
            %       using data from SSTV <superSpaceTimeVolumeIn> in frames defined by <frameGroup>, which is a vector of frame indices (e.g. 1:15).
            
            if (nargin<1)
                CorrespondingPointListInstance.data = [];
            elseif (nargin == 1)
                
                if ( strcmpi(class(varargin{1}) ,'CorrespondingPointList') ) %Initialize from another CPL
                    
                    nCamGroups = size(varargin{1},3);
                    
                    CorrespondingPointListInstance(nCamGroups) = CorrespondingPointList();
                    
                    for iterCamGroup = 1:nCamGroups
                        CorrespondingPointListInstance(iterCamGroup).data = varargin{1}(iterCamGroup).data;
                        CorrespondingPointListInstance(iterCamGroup).normalizeCoordinates = varargin{1}(iterCamGroup).normalizeCoordinates;
                        CorrespondingPointListInstance(iterCamGroup).calcFMethod = varargin{1}(iterCamGroup).calcFMethod;
                    end
                    
                else
                    disp('Wrong Class at input of CPL Constructor');
                    keyboard;
                end
                
            elseif (nargin == 3) %Initialize from Tracks, SSTV and Frame Limits
                if ( strcmpi(class(varargin{1}) ,'SuperTrackList') && strcmpi(class(varargin{2}) ,'SuperSpaceTimeVolume') )
                    STLIn = varargin{1};
                    SSTVIn = varargin{2};
                    frameGroup = varargin{3};
                    
                    nCamGroups = size(STLIn.data,1);
                    
                    if (nCamGroups ~= size(SSTVIn, 4)-1) %** 1 For Groups = Pairs**%
                        disp('Contradicting Number of Camera Groups between STL and SSTV!');
                        keyboard;
                    end
                    
                    CorrespondingPointListInstance(nCamGroups) = CorrespondingPointList();
                    
                    for iterCamGroup = 1:nCamGroups
                        
                        camGroup = iterCamGroup : (iterCamGroup+1); %**Passing two views for corresponding SuperTrack**%
                        
                        CorrespondingPointListInstance(iterCamGroup).data = ...
                            CorrespondingPointListInstance.ExtractPointsFromTracks(STLIn.data(iterCamGroup,:), SSTVIn(:, :, frameGroup, camGroup));
                    end
                    
                else
                    disp('Wrong Class at input of CPL Constructor');
                    keyboard;
                end
                
            else
                disp('Wrong Number of parameters at input of CPL Constructor');
                keyboard;
                
            end %nargin check
        end %constructor
        
        %METHODS
        
        function fundMat = CalcFundMat(self, varargin)
            %   FUNDAMENTAL MATRIX CALCULATION
            %
            %       fundMat = CPL.CalcFundMat() - Calculate the Fundamental Matrix from the List of Matching Points using the method set in <calcFMethod>.
            %
            %       fundMat = CPL.CalcFundMat(method) - Calculate the Fundamental Matrix from the List of Matching Points, using the method <method>,
            %           which can be 'LS' for simple Least Squares, 'RANSAC', or 'IRLS' for Iteratively Reweighted Least Squares (default).
            %
            %       fundMat = CPL.CalcFundMat(method, normalize) - Same as above, <normalize> is a binary variable stating whether to perform coordinate normalization
            %           or not. (default = 1)
            %
            %       In case of more than two cameras (CPL is an array, and the 4th dimention of the SSTV is of the appropriate size - for pairs, one camera more than the number of
            %           groups), fundMat is a 3D array, having a 3x3 layer for each camera group.
            
            if (nargin >= 2) %Method given at input
                
                if ( strcmpi(varargin{1} ,'IRLS') || strcmpi(varargin{1} ,'RANSAC') || strcmpi(varargin{1} ,'LS') )
                        self.calcFMethod = varargin{1};
                    
                else
                    disp('Wrong Method at input of CalcFundMat');
                    keyboard;
                end
            end
            
            if (nargin == 3)
                if ( (varargin{2} == 0) || (varargin{2} == 1) )
                    self.normalizeCoordinates = varargin{2};
                else
                    disp('Wrong Coordinate Normalization Flag at input of CalcF (only 0 or 1 allowed)');
                    keyboard;
                end
            end
            
            if (nargin > 3)
                disp('Wrong number of parameters at input of CalcF (see help for allowed parameters)');
                keyboard;
            end
            
            nCamGroups = size(self,1);
            
            fundMat = zeros(3,3,nCamGroups);
            
            for iterCamGroup = 1:nCamGroups
                %****CURRENTLY ONLY FOR PAIRS*****%
                
                if strcmpi(self.calcFMethod ,'IRLS')
                    fundMat = self.IRLSF(self(iterCamGroup,:).data(:,1:2), self(iterCamGroup,:).data(:,3:4));
                    
                elseif strcmpi(self.calcFMethod ,'RANSAC')
                    [fundMat(:,:,iterCamGroup), ~, ~, ~] = self.ransacF(self(iterCamGroup,:).data(:,1:2), self(iterCamGroup,:).data(:,3:4));
                    
                elseif strcmpi(self.calcFMethod ,'LS')
                    fundMat(:,:,iterCamGroup) = self.calcF(self(iterCamGroup,:).data(:,1:2), self(iterCamGroup,:).data(:,3:4));
                end
                    
            end
            
        end
        
        function self = AddPointsFromAnotherSlice(self, superTrackListIn, superSpaceTimeVolumeIn, frameGroup)
            %   CPL.AddPointsFromAnotherSlice(superTrackListIn, superSpaceTimeVolumeIn, frameGroup)
            %   	Add points extracted from STL <superTrackListIn> (assuming STL is already matched),
            %       using data from SSTV <superSpaceTimeVolumeIn> in frames defined by <frameGroup>, which is a vector of frame indices (e.g. [1, 15])
            
            if ( ~strcmpi(class(superTrackListIn) ,'SuperTrackList') || ~strcmpi(class(superSpaceTimeVolumeIn) ,'SuperSpaceTimeVolume') )
                disp('Wrong Class at input of CPL point addition');
                keyboard;
            else
                
                nCamGroups = size(superTrackListIn.data,1);
                
                if (nCamGroups ~= size(superSpaceTimeVolumeIn, 4)-1) %** 1 For Groups = Pairs**%
                    disp('Contradicting Number of Camera Groups between STL and SSTV!');
                    keyboard;
                end
                
                if ( size(self, 1) ~= nCamGroups)
                    disp('Contradicting Number of Camera Groups between old CPL and new SSTV!');
                    keyboard;
                end
                
                for iterCamGroup = 1:nCamGroups
                    
                    camGroup = iterCamGroup : (iterCamGroup+1); %**Passing two views for corresponding SuperTrack**%
                    
                    self(iterCamGroup).data = [self(iterCamGroup).data; ...
                        self.ExtractPointsFromTracks(superTrackListIn.data(iterCamGroup,:), superSpaceTimeVolumeIn(:, :, frameGroup, camGroup))];
                end
                
            
            end
            
            
        end
        
        function self = RemoveDoublePoints(self)
            %CPL.RemoveDoublePoints() - Remove points occuring more than once in the list
            
            for iterCamGroup = 1:length(self)
                [~, noDuplicates, ~] = unique(self(iterCamGroup).data,'rows', 'stable');
                self(iterCamGroup).data = self(iterCamGroup).data(noDuplicates,:);
            end
            
        end
        
    end %Public Methods
    
%     methods (Static)
%         function distances = CalcGeometricDistances(fundMat, P1, P2)
%             %geometricDistances = CPL.CalcGeometricDistances(fundMat, P1, P2)
%             %
%             %Calculates the geometric desitances between epipolar lines
%             %   calculated from fundMat and P1, and the corresponding points
%             %   P2. Normalization must be performed outside of this
%             %   function.
%             
%             fp = fundMat * P1';
%             distances = abs(diag( P2*(fp./(repmat(abs(sqrt(fp(1,:).^2+fp(2,:).^2)),3,1)+eps))) );
%             
%         end
%     end
    
    methods (Access = private)
        
        function points = ExtractPointsFromTracks(self, STLRow, SSTVSlice)
            %points = ExtractPointsFromTracks(STLLayer, SSTV2FrameSlice)
            %
            %Extracts Interest Points from List of Tracks in STLLayer (single layer of STL array), using full data from the Slice of SSTV <SSTVSlice>.
            
            sliceSize = size(SSTVSlice);
            
            nTracks = size(STLRow{1},1);
            
            nFrames = sliceSize(3);
            
            nCamsInGroup = sliceSize(4);
            
            points = zeros(nTracks*nFrames, 2*sliceSize(4));
            
            for iterCamInGroup = 1:nCamsInGroup
                
                curCamGroup = (iterCamInGroup-1)*nCamsInGroup + (1:nCamsInGroup);
                
                for iterTrack = 1:nTracks  %for each track
                    
                    curTrack = STLRow{iterCamInGroup};
                    
                    curTrackGroupBeginning = (iterTrack-1)*nFrames;
                    
                    %create mask of track location
                    curTrackMask = zeros(sliceSize(1:2));
                    curTrackMask(...
                        round(curTrack(iterTrack).BoundingBox(2)) + (0:curTrack(iterTrack).BoundingBox(4)),...
                        round(curTrack(iterTrack).BoundingBox(1)) + (0:curTrack(iterTrack).BoundingBox(3))) = 1;
                    
                    curTrackMask = repmat(curTrackMask,[1 1 nFrames]);
                    
                    SingleViewSlice = SSTVSlice(:,:,:,iterCamInGroup).*curTrackMask;
                    
                    
                    for iterFrame = 1:nFrames %for each point extracted from track
                        
                        curPoint = curTrackGroupBeginning + iterFrame;
                        
                        
                        [x, y] = self.extractMaximaLocation(SingleViewSlice(:,:,iterFrame));
                        
                        points( curPoint, curCamGroup) = [x y];
                        
                        
                    end
                    
                end
                
            end
            
        end
        
        function [x, y] = extractMaximaLocation(~, imgIn)
            %[y x] = extractMaximaLocation(imgIn)
            %
            %Extract the location of the (first) maximal value in the image <imgIn>
            
            [y, x] = find(imgIn == max(imgIn(:)),1);
            
        end
        
        function F = calcF(self, p1, p2, weightsOfPoints)
            %F = calcF(p1, p2, weightsOfPoints)
            %8 pt. alg. Coordinate normalization dependant on the field in self
            %
            %p1, p2 - nx2 vectors of corresponding point coordinates
            %weightsOfPoints (optional) - the weights of the points in vector form (nPointsx1). default = ones.
            
            nPoints = size(p1,1);
            
            if ~exist('weightsOfPoints', 'var')
                weightsOfPoints = ones(nPoints,1);
            end
            
            %%%Coordinates Normalization%%%
            if self.normalizeCoordinates
                %mean point
                meanP1 = mean(p1,1);
                meanP2 = mean(p2,1);
                
                %stDev of points
                stdP1 = mean(sqrt(sum((p1 - repmat(meanP1,nPoints,1)).^2,2)));
                stdP2 = mean(sqrt(sum((p2 - repmat(meanP2,nPoints,1)).^2,2)));
                
                %Normalization Matrix
                T1 = [ [eye(2), -meanP1']/stdP1; [0 0 1] ];
                T2 = [ [eye(2), -meanP2']/stdP2; [0 0 1] ];
                
                normP1 = (T1*[p1, ones(nPoints, 1)]')';
                normP2 = (T2*[p2, ones(nPoints, 1)]')';
                
            else
                normP1 = p1;
                normP2 = p2;
            end
            
            
            %%%8 Point Alg.%%%
            A = [normP1(:,1).*normP2(:,1), normP1(:,1).*normP2(:,2), normP1(:,1), normP1(:,2).*normP2(:,1),...
                normP1(:,2).*normP2(:,2), normP1(:,2), normP2(:,1), normP2(:,2), ones(nPoints,1)];
            
            A = diag(weightsOfPoints)*A;
            
            AtA = A'*A;
            
            [~, ~, V] = svd(AtA);
            
            F = reshape(V(:,end), 3, 3);
            
            [U, S, V] = svd(F);
            S(end) = 0;
            
            F = U*S*V';
            
            if self.normalizeCoordinates
                %%%Return to the original Coordinates
                F = T2'*F*T1;
                F = F/norm(F,2);
            else
                F = F/norm(F,2);
            end
            
        end
        
        function [bestF, nInliers, medianError, inliersOut] = ransacF(self, p1,p2,f,s, errorThreshold, goal)
            % CPL.ransacF(p1, p2, f, s, errorThreshold, goal)
            %
            % A function to calculate the fundumental matrix F.
            % The function uses RANSAC for calculation of F.
            % Arguments:
            %	p1,p2                   Nx2 arrays that contain pixel coordinates of N  matching points.
            %	f                       the fraction of inliers in the data (default = 0.6).
            %	s                       the desired probability of success (default = 0.99).
            %   errorThreshold          the threshold fr defining the point as an outliers (in pixels from the epipolar line). (default = 3)
            %   goal                    'inliers' for maximal number of inliers (default), or 'error' for minimal median error
            %
            % Optional Outputs:
            %   nInliers - the number of inliers found by the method
            %   medianError - the median error of the inliers
            %   inliersOut - the nInliersx4 array of the inlier points
            
            %%
            % calculate number of retries
            if ~exist('f', 'var')
                f = 0.6;	% fraction of inliers
            end
            
            if ~exist('s', 'var')
                s = 0.99;	% probability of success
            end
            
            if ~exist('errorThreshold', 'var')
                errorThreshold = 3;	% acceptable error in pixels from epipolar line
            end
            
            if ~exist('goal', 'var')
                goal = 'inliers';	% acceptable error in pixels from epipolar line
            end
            
            if ~(strcmpi(goal, 'error')||strcmpi(goal, 'inliers'))
                disp('Wrong goal!')
                keyboard;
            end
            
            nPairs = 8;
            k = ceil ( log(1-s) / log ( 1 - f^nPairs ) );
            
            % number of matching points
            nPoints = size (p1,1);
            
            %%%Coordinates Normalization%%%
            if self.normalizeCoordinates
                %mean point
                meanP1 = mean(p1,1);
                meanP2 = mean(p2,1);
                
                %stDev of points
                stdP1 = mean(sqrt(sum((p1 - repmat(meanP1,nPoints,1)).^2,2)));
                stdP2 = mean(sqrt(sum((p2 - repmat(meanP2,nPoints,1)).^2,2)));
                
                %Normalization Matrix
                T1 = [ [eye(2), -meanP1']/stdP1; [0 0 1] ];
                T2 = [ [eye(2), -meanP2']/stdP2; [0 0 1] ];
                
                normP1 = (T1*[p1, ones(nPoints, 1)]')';
                normP2 = (T2*[p2, ones(nPoints, 1)]')';
                
                normP1 = normP1(:,1:2);
                normP2 = normP1(:,1:2);
                
            else
                normP1 = p1;
                normP2 = p2;
            end
            
            %%
            bestScore = Inf;
            maxNInliers = 0;
            %do k times
            for iterRANSAC=1:k
                
                possibleInliers = randperm(nPoints);
                possibleF = self.calcF(normP1(possibleInliers(1:nPairs),:) , normP2(possibleInliers(1:nPairs),:));
                
                fp = possibleF * [normP1'; ones(1,nPoints)];
                allErrors = abs( diag( [normP2, ones(nPoints,1)]*fp/abs(sqrt(fp(1)^2+fp(2)^2)) ) );
                currentInliers  = (allErrors<errorThreshold);
                
                if sum(currentInliers) > (nPairs-1)
                    currentF = self.calcF(normP1(currentInliers,:) , normP2(currentInliers,:));
                    allCurrentErrors = diag( abs( [normP2, ones(nPoints,1)] * currentF * [normP1'; ones(1,nPoints)] ) );
                    
                    score = median( allCurrentErrors(currentInliers) );
                    
                    if strcmpi(goal, 'error')
                        condition = (score < bestScore);
                    else
                        condition = (sum(currentInliers) > maxNInliers);
                    end
                    
                    if condition
                        nInliers = sum (currentInliers);
                        bestScore = score;
                        maxNInliers = sum (currentInliers);
                        
                        if self.normalizeCoordinates
                            %%%Return to the original Coordinates
                            bestF = T2'*currentF*T1;
                            bestF = bestF/norm(bestF,2);
                        else
                            bestF = currentF;
                            bestF = bestF/norm(bestF,2);
                        end
                        
                        medianError = score;
                        inliersOut = [p1(currentInliers,:) p2(currentInliers,:)];
                    end
                    
                end
                
            end
        end
        
        function fundMat = IRLSF(self, p1,p2, numIterations)
            % CPL.IRSLF(p1, p2, numIterations)
            %
            % A function to calculate the fundumental matrix F.
            % The function uses IRLS (Iteratively Reweighted Least Squares for calculation of F.
            % Arguments:
            %	p1,p2 are Nx2 arrays that contain pixel coordinates of N  matching points.
            %   numIterations is the maximal number of iterations the IRLS will run (default = 20)
            
            if ~exist('numIterations', 'var')
                numIterations = 20;
            end
            
            % number of matching points
            nPoints = size (p1,1);
            
            if self.normalizeCoordinates
                %mean point
                meanP1 = mean(p1,1);
                meanP2 = mean(p2,1);
                
                %stDev of points
                stdP1 = mean(sqrt(sum((p1 - repmat(meanP1,nPoints,1)).^2,2)));
                stdP2 = mean(sqrt(sum((p2 - repmat(meanP2,nPoints,1)).^2,2)));
                
                %Normalization Matrix
                T1 = [ [eye(2), -meanP1']/stdP1; [0 0 1] ];
                T2 = [ [eye(2), -meanP2']/stdP2; [0 0 1] ];
                
                normP1 = (T1*[p1, ones(nPoints, 1)]')';
                normP2 = (T2*[p2, ones(nPoints, 1)]')';
                
            else
                normP1 = [p1; ones(nPoints, 1)];
                normP2 = [p2; ones(nPoints, 1)];
            end
            
            %IRLS
            display('Starting IRLS...');
            
            self.pointWeights = ones(nPoints,1);
            
            for iterIRLS=1:numIterations
                
                fundMat = self.calcF(normP1(:,1:2) , normP2(:,1:2) , self.pointWeights);
                
                geometricDistances = CalcGeometricDistances(fundMat, normP1, normP2);
                
                self.pointWeights = min(1./geometricDistances, 1./median(geometricDistances));
            end
            
            display('IRLS Done.');
            
            if self.normalizeCoordinates
                %%%Return to the original Coordinates
                fundMat = T2'*fundMat*T1;
            end
            
            fundMat = fundMat/norm(fundMat,2);
            
        end
        
    end % private methods
    
end % class