classdef StereoTrackList < handle
    %StereoTrackList is a class defining the List holding the Tracks detected in the StereoSlice, for the PTV Point Matching Algorithm
    %
    %Data Structure: A Cell Array representing cameraView pairs (see below) along the first direction, and the Slices of the Pair in the second direction.
    %                In each cell, is an 1xN array of structs, each struct representing the parameters of one of the N tracks detected by Connecteed Compunents
    %                Labeling.
    %
    %Properties: cannyThreshold = 0.3)
    %            AngleLimit (default = 90)
    %            ImageLength = (default: set from length of slice image)
    %            PreferenceOfPositionOverOrientation (default = 3)
    %
    %CONSTRUCTORS:
    %   stereoTrackListInstance = SteroTrackList() - Create empty object of StreoTrackList Class
    %
    %   stereoTrackListInstance = SteroTrackList(sTLInstanceOld) - Initialize from another STL Instance <sTLInstanceOld>
    %
    %   stereoTrackListInstance = SteroTrackList(superSliceIn) - Initialize by detecting tracks in StereoSlice <steroSliceIn>
    %
    %METHODS:
    %   sTLInstance = STL.Match() - Match tracks over the current Stereo List and return the matching tracks in a new list <sTLInstance>
    %       **In case of several cameras - matching is between all consecutive pairs (1-2, 2-3, 3-4, ...).
    %           Then, STLInstance is in itself an array containing in each place the Track List of the corresponding consecutive pair.
    %           **Pay attention to inverting images where necessary...
    %
    % % % %     %   track = STL.GetTrack(trackNum, cameraNum, direction) - Retrieve track number <trackNum> in <cameraNum> camera view.
    % % % %     %       In case of more than two Camera Views, <direction> is used to describe the camera position in the pair. <direction> can be 'forward' (default) or
    % % % %     %       'backward'.
    %
    %Simon Kolotov, Ver 2.0, Spring 2013
    
    properties (SetAccess = public, GetAccess = public)
        data
        
        cannyThreshold = 0.3;
        
        %Normalization Parameters for the Matching Method
        AngleLimit = 90;
        ImageLength = [];
        PreferenceOfPositionOverOrientation = 3;
    end %Public Properties
    
    methods
        %CONSTRUCTORS
        function stereoTrackListInstance = StereoTrackList(varargin)
            %   stereoTrackListInstance = StereoTrackList() - Create empty object of StreoTrackList Class
            %
            %   stereoTrackListInstance = StereoTrackList(sTLInstanceOld) - Initialize from another STL Instance <sTLInstanceOld>
            %
            %   stereoTrackListInstance = StereoTrackList(stereoSliceIn) - Initialize by detecting tracks in StereoSlice <steroSliceIn>, using
            %       default canny Threshold
            %
            %   stereoTrackListInstance = StereoTrackList(stereoSliceIn, cannyThresholdIn) - Initialize by detecting tracks in StereoSlice <steroSliceIn>, using
            %       default canny Threshold
            %
                        
            if (nargin<1)
                stereoTrackListInstance.data = [];
                return
            elseif (nargin == 1)
                
                if ( isequal(class(varargin{1}) ,'StereoTrackList') ) %Initialize from another STL
                    
                    nCamPairs = size(varargin{1},3);
                    
                    stereoTrackListInstance(nCamPairs) = StereoTrackList();
                    for iterCamPairs = 1:nCamPairs
                        stereoTrackListInstance.data = varargin{1}.data;
                        stereoTrackListInstance.cannyThreshold = varargin{1}.cannyThreshold;
                        stereoTrackListInstance.AngleLimit = varargin{1}.AngleLimit;
                        stereoTrackListInstance.ImageLength = varargin{1}.ImageLength;
                        stereoTrackListInstance.PreferenceOfPositionOverOrientation = varargin{1}.PreferenceOfPositionOverOrientation;
                        return;
                    end
                    
                elseif ( strcmpi(class(varargin{1}) ,'SuperSlice') ) %Initialize by Track Detection With Default Canny Threshold
                    sliceIn = varargin{1};
                    
                    stereoTrackListInstance.ImageLength = max(size(sliceIn)); %Length of the image {assumption - nCameraViews << length(image) }
            
                    for iterCamView = 1:size(sliceIn,3)
                        stereoTrackListInstance.data{1,iterCamView} = stereoTrackListInstance.detectTracksInSlice(sliceIn(:,:,iterCamView));
                    end
                    
                    
                else
                    disp('Wrong Class at input of STL Constructor');
                    keyboard;
                end
                
            elseif (nargin ==2) %Initialize by Track Detection With given Canny Threshold
                if ( strcmpi(class(varargin{1}) ,'SuperSlice') ) 
                    sliceIn = varargin{1};
                    
                    stereoTrackListInstance.ImageLength = max(size(sliceIn)); %Length of the image {assumption - nCameraViews << length(image) }
                    
                    stereoTrackListInstance.cannyThreshold = varargin{2};
            
                    for iterCamView = 1:size(sliceIn,3)
                        stereoTrackListInstance.data{1,iterCamView} = stereoTrackListInstance.detectTracksInSlice(sliceIn(:,:,iterCamView));
                    end
                    
                    
                else
                    disp('Wrong Class at input of STL Constructor');
                    keyboard;
                end
            else
                disp('Too many parameters at input of STL Constructor');
                keyboard;
                
            end %nargin check

        end %constructor
        
        %METHODS
        
        function STLInstance = Match(self)
            %sTLInstance = STL.Match()
            %
            %Match tracks over the current Stereo List and return the matching tracks in a new list <sTLInstance>
            %
            %**In case of several cameras - matching is between all consecutive pairs (1-2, 2-3, 3-4, ...).
            %    Then, STLInstance is in itself an array containing in each place the Track List of the corresponding consecutive pair.
            %    **Pay attention to inverting images where necessary...
            
            if ( size(self.data,1) > 1 ) %Tracks already parsed ###Note - partial check, only for more than 2 camera views
                disp('Tracks already matched (only arrays of 1xN are acceptable)');
                keyboard;
            end
            
            nCorrespondingPairs = size(self.data,2)-1;
            nCameras = size(self.data,2);
            
            STLInstance = StereoTrackList(); %Initialize Array (for more than two camera views)
            STLInstance.data = cell(nCorrespondingPairs, nCameras);
            
            for iterCamViews = 1:nCorrespondingPairs
                
                %Correspond forward
                pairsForward = self.findSimilarDescriptors(self.data{iterCamViews}, self.data{iterCamViews + 1});
                %Correspond backward
                pairsBackward = self.findSimilarDescriptors(self.data{iterCamViews+1}, self.data{iterCamViews});
                
                %Find Double Correspondances
                doublePairsN = [];
                doublePairsNPlus1 = [];
                for iterPairs = 1:size(pairsForward,1)

                    curPair = repmat(fliplr(pairsForward(iterPairs,:)), size(pairsBackward,1), 1); %
                    locDouble = find( sum( ~(pairsBackward-curPair) , 2) == 2); %locate those pairs that appear in both lists
                    
                    doublePairsN = [doublePairsN; pairsBackward(locDouble, 2)]; %#ok (double pairs change size every iteration)
                    doublePairsNPlus1 = [doublePairsNPlus1; pairsBackward(locDouble, 1)]; %#ok (double pairs change size every iteration)
                end
                
                STLInstance.data{iterCamViews, 1} = self.data{iterCamViews}(doublePairsN);
                STLInstance.data{iterCamViews, 2} = self.data{iterCamViews+1}(doublePairsNPlus1);
                
            end
            
        end
        
% % %         function track = GetTrack(self, trackNum, cameraNum, varargin)
% % %             %track = STL.GetTrack(trackNum, cameraNum)
% % %             %
% % %             %Retrieve track number <trackNum> in <cameraNum> camera view.
% % %             %    In case of more than two Camera Views, <direction> is used to describe the camera position in the pair.
% % %             %    <direction> can be 'forward' (default) or 'backward'.
% % %             nCorrespondingPairs = size(self, 2);
% % %             
% % %             if (nCorrespondingPairs < 1)
% % %                 disp('Cannot extract track, TrackList is empty!');
% % %                 keyboard;
% % %             
% % %             elseif (nCorrespondingPairs == 1)   %one pair
% % %                 cameraPair = 1;
% % %             
% % %             else    %more than one pair
% % %                 if (nargin < 4) %direction is not specified
% % %                     
% % %                     direction = 'forward'; %Default - Forward
% % %                     
% % %                 elseif (nargin ==4)
% % %                     
% % %                     direction = varargin{1};
% % %                     
% % %                     if strcmpi(direction, 'forward')
% % %                         cameraPair = cameraNum;
% % %                         cameraNum = 1;
% % %                     elseif strcmpi(direction, 'forward')
% % %                         cameraPair = cameraNum - 1;
% % %                         cameraNum = 2;
% % %                     else
% % %                         disp ('Unknown Direction');
% % %                         keyboard;
% % %                     end
% % %                     
% % %                 else
% % %                     disp ('Too Many Input Variables');
% % %                     keyboard;
% % %                 end
% % %             end
% % %             track = self(cameraPair).data{cameraNum}(trackNum);
% % %         end
        
        
    end %Public Methods
    
    methods (Access = private)

        function trackList = detectTracksInSlice(self, sliceImage)
            %trackList = STL.detectTracksInSlice(sliceImage)
            %
            %performs track detection in <sliceImage> based on canny edge detection and connected components labeling.
            
            imgEdge = edge(sliceImage, 'canny', self.cannyThreshold);
            
            trackList = regionprops(imgEdge, sliceImage, 'all');
            
        end
        
        function pairs = findSimilarDescriptors(self, tracks, refTracks) %**Only for Matching paris of trackLists**%
            %pairs = STL.findSimilarDescriptors(regionsOfInterest, regionsOfInterestRef)
            
            normalizationOfParameters = [ [[1 1] self.PreferenceOfPositionOverOrientation*[1 1]]/self.ImageLength 1/self.AngleLimit];
            %Major axis, Minor axis, centroidX, centroidY, Orientation
                        
            
            %Create Reference Descriptors
            refCentroids = [refTracks(:).Centroid]';
            refCentroidsX = refCentroids(1:2:end);
            refCentroidsY = refCentroids(2:2:end);
            
            refDescriptor = [[refTracks(:).MajorAxisLength]' [refTracks(:).MinorAxisLength]'...
                refCentroidsX refCentroidsY [refTracks(:).Orientation]'].*repmat(normalizationOfParameters,size(refTracks,1),1);
            

            nTracks = length(tracks);
            leavingLocs = [];
            closestTrack = inf(nTracks,2);
            
            for iterTrack = 1:nTracks
                %Create Copies of Current Descriptor
                curCentroid = [tracks(iterTrack).Centroid]';
                curCentroidX = curCentroid(1);
                curCentroidY = curCentroid(2);
                
                curDescriptor = [tracks(iterTrack).MajorAxisLength tracks(iterTrack).MinorAxisLength...
                    curCentroidX curCentroidY tracks(iterTrack).Orientation].*normalizationOfParameters;
                curDescriptor = repmat(curDescriptor,size(refDescriptor,1), 1);
                
                %Calculate Euclidian Distances between current descriptor and reference descriptors
                curEuclidianDist = sum((refDescriptor-curDescriptor).^2,2);
                
                %Find Closest Descriptors to every descriptor
                minDistLoc = find(curEuclidianDist==min(curEuclidianDist),1);
                closestTrack(iterTrack, :) = [minDistLoc min(curEuclidianDist)];
                
                redun = find(closestTrack(:,1)==minDistLoc); %Other Tracks have this neighbor
                leavingLocs = [leavingLocs; redun(closestTrack(redun,2) ~= min(closestTrack(redun,2)))]; %#ok (LeavingLocs grows every iteration)
            end
            
            stayingLocs = (1:length(closestTrack))';
            stayingLocs(leavingLocs) = []; %Remove Tracks That have no neighbor of their own
            
            closestTrack(leavingLocs,:) = []; %Remove Corresponding Neighbors
            
            pairs = [stayingLocs closestTrack(:,1)]; %return pairs of tracks having their own neighbors and those neighbors
            
        end
   end % private methods

end % class