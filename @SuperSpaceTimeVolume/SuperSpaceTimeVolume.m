classdef SuperSpaceTimeVolume < handle
    %SuperSpaceTimeVolume is a class defining the super space-time volume for the PTV Point Matching Algorithm
    %
    %DATA ACCESS:
    %   SSTV(X, Y, Frame, CamView), where:
    %                   (X,Y) are the coordiantes in a single image,
    %                   Frame is the frame number,
    %                   CamView is the Index of the camera in the Camera Array
    %
    %FundMat: 3x3xNCamGroups matrix, each layer representing the fundamental matrix of each group of cameras in the array.
    %   This variable is not initialized internally, but calculated from a corresponding point list.
    %
    %K, R, t are the Internal Parameters Matrix (3x3), and Rotational (3x3) and translational (3x1)
    %   components of the External Parameters Matrix.
    %
    %CONSTRUCTORS:
    %   superSpaceTimeVolumeInstance = SuperSpaceTimeVolume() - Create empty object of SuperSpaceTimeVolume Class
    %
    %   SpaceTimeVolumeInstance = SuperSpaceTimeVolume(SSTVInstanceOld) - Initialize STV from another STV Instance <stvInstanceOld>
    %
    %   superSpaceTimeVolumeInstance = SuperSpaceTimeVolume(pathIn, nameDescriptors) - Initialize SSTV from path <pathIn> for image
    %   names beginning with {<nameDescriptor1>; <nameDescriptor2>; ...} (e.g. {'Left_1'; 'Right_1'})
    %
    %METHODS:
    %   SSTV.GetImage(frameNum, cameraNum) - retrieve a specific image at frame <frameNum> from <cameraNum> view
    %
    %   F = SSTV.ResampleF(resampleRate) - Upsample or downsample the fundamental matrix of SSTV by a factor <resampleRate>, and return it.
    %
    %   SSTV_D = SSTV.DownSample(downSampleRate) - Create a new SSTV by downsampling every frame of the SSTV by <downsampleRate>
    %
    %   SSTV.EpipolarDistances(frameNumber, groundTruth) - Displays frame <frameNumber> of the sequence from all views,
    %       and allows to check distances from the epipolar lines of selected points to points selected in the corresponding image.
    %       Optional paramter groundTruth allows to overlay ground truth points (given in a Nx4 matrix) and the median error of
    %       the fundamental matrix over those points.
    %
    %   SSTV.RelativeCalibration(method, points1, points2, f0, Cx0, Cy0) - Performs relative calibration according to the
    %       PTV Point Matching algorithm, setting the K, R and t parameters of the SSTV.
    %       <method> is 'full' for full grid search, 'sparse2full' for sparse to full search, or 'separate' for separate grid search.
    %       <points1> and <points2> are two lists of points (each a Nx2 array of pixel coordinates) to be used for optimization, 1 being for the left camera,
    %       and 2 for the right. **CURRENTLY ONLY FOR 2 CAMERAS**
    %       <f0>, <Cx0>, and <Cy0> are the starting parameters for the search.
    %
    %   medianDistance = SSTV.CalculateEpipolarDistances(groundTruth)
    %       Similar to EpipolarDistances, without the display. Instead of
    %       displaying the results graphically, only the median error of the fundamental matrix over the GroundTruth points is calculated and returned.
    %
    %Simon Kolotov, Ver 2.0, Spring 2013
    
    properties (SetAccess = public, GetAccess = public)
        data
        
        fundMat
        
        K, R, t
        
        colors = 'grcmywk';
        shapes = 'osdv^<>ph';
        
    end %Public Properties
    
    methods
        %CONSTRUCTORS
        function SuperSpaceTimeVolumeInstance = SuperSpaceTimeVolume(varargin)
            %SpaceTimeVolumeInstance = SpaceTimeVolume() - Create empty object of SpaceTimeVolume Class
            %
            %SpaceTimeVolumeInstance = SpaceTimeVolume(STVInstanceOld) - Initialize STV from another STV Instance <stvInstanceOld>
            %
            %SpaceTimeVolumeInstance = SpaceTimeVolume(pathIn, nameDescriptors) - Initialize SSTV from path <pathIn>
            %   (ending with "<Sequence Name>\") for image names beginning with {<nameDescriptor1>; <nameDescriptor2>; ...}  (e.g. {'Left_1'; 'right1'})
            
            if (nargin== 1)
                
                if (strcmpi(class(varargin{1}) ,'SuperSpaceTimeVolume')) %Initialize from another SSTV
                    SuperSpaceTimeVolumeInstance.data = varargin{1}.data;
                    SuperSpaceTimeVolumeInstance.fundMat = varargin{1}.fundMat;
                else
                    disp('Wrong Class at input of SSTV Constructor');
                    keyboard;
                end
                
            elseif (nargin == 2)
                
                if ~exist(varargin{1}, 'dir')
                    disp('Wrong Path at input of SSTV Constructor');
                    keyboard;
                end
                
                pathIn = varargin{1};
                
                % -_-_-_-_Locate the images-_-_-_-_%
                nameDescriptors = varargin{2};
                
                nViews = size(nameDescriptors,1);
                
                imageNames = dir([pathIn  nameDescriptors{1} '*.tif']);
                nImages = length(imageNames); %Number of such images in the directory
                
                if ~nImages
                    disp(['No Images beginning with ' nameDescriptors{1} ' Found at ' pathIn '!']);
                    keyboard;
                end
                
                imgSize = size(imread([pathIn imageNames(1).name]));
                
                SuperSpaceTimeVolumeInstance.data = NaN([imgSize(1:2) nImages nViews]);
                
                for iterCamViews = 1:nViews
                    
                    imageNames = dir([pathIn  nameDescriptors{iterCamViews} '*.tif']);
                    nImages = length(imageNames); %Number of such images in the directory
                    %!!!!ONLY LOAD, NO CONSISTENCY CHECK!!!!
                    
                    if ~nImages
                        disp(['No Images beginning with ' nameDescriptors{iterCamViews} ' Found at ' pathIn '!']);
                        keyboard;
                    end
                    
                    
                    
                    % -_-_-_-_Start Aggregation-_-_-_-_%
                    disp(['Starting image aggregation on ' nameDescriptors{iterCamViews} '...']);
                    for iterImg = 1:nImages  %Iterate through images
                        %% Load Images
                        imgIn = double(imread([pathIn imageNames(iterImg).name]));
                        
                        %% Aggregate Images
                        SuperSpaceTimeVolumeInstance.data(:,:,iterImg,iterCamViews) = imgIn(:,:,1);
                        %!!!!GRAYSCALE IMAGES ONLY. NO DOWNSAMPLE, NOT QUANTIZATION!!!!
                        
                        %I'm Alive
                        if ~mod(iterImg, 75)
                            disp(['image ' num2str(iterImg) ' of ' num2str(nImages) ' Done.']);
                        end
                    end
                    
                end %iterCamViews
                
                disp('Image aggregation Done.');
                
            elseif (nargin > 2)
                disp('Too many parameters at input of STL Constructor');
                keyboard;
                
            end %nargin check
            
        end %constructor
        
        %METHODS
        function img = GetImage(self, frameNumber, cameraNumber)
            %SSTV.GetImage(frameNum, cameraNum)
            %
            %retrieve a specific image at frame <frameNum> from <cameraNum> view
            
            if (frameNumber > size(self.data,3))
                disp(['Frame ' num2str(frameNumber) ' out of bounds of the SSTV']);
                keyboard;
            elseif (cameraNumber > size(self.data,4))
                disp(['Camera ' num2str(frameNumber) ' out of bounds of the Camera Array']);
                keyboard;
            else
                img = self.data(:,:,frameNumber, cameraNumber);
            end
        end
        
        function superSpaceTimeVolumeInstance = Downsample(self, downsampleRate)
            %SSTV_D = SSTV.DownSample(downSampleRate)
            %
            %Create a new SSTV by downsampling every frame of the SSTV by <downsampleRate>
            
            superSpaceTimeVolumeInstance = SuperSpaceTimeVolume(self);
            
            superSpaceTimeVolumeInstance.data = [];
            
            if ~isempty(self.fundMat)
                superSpaceTimeVolumeInstance.fundMat = self.ResampleF(downsampleRate);
            end
            
            for iterCamViews = 1:size(self.data,4)
                disp(['Starting Downsampling on View ' num2str(iterCamViews) '...']);
                
                superSpaceTimeVolumeInstance.data(:,:,:,iterCamViews) = imresize(self.data(:,:,:,iterCamViews),1/downsampleRate);
            end
            disp('Done.');
            
        end
        
        function resF = ResampleF(self, resampleRate)
            %F = SSTV.ResampleF(resampleRate)
            %
            %Upsample or downsample the fundamental matrix of SSTV by a factor <resampleRate>, and return it.
            
            resamplingMatrix = [1 0 0 ; 0 1 0; 0 0 1/resampleRate]*resampleRate;

            resF = resamplingMatrix*self.fundMat*resamplingMatrix;
        end
        
        function EpipolarDistances(self, frameNumber, groundTruth)  %CURRENTLY FOR PAIRS ONLY!!
            % SSTV.EpipolarDistances(frameNumber, groundTruth)
            % Displays frame <frameNumber> of the sequence from all views, and allows to check distances from the epipolar lines of
            %       selected points to points selected in the corresponding image.
            % Optional paramter groundTruth allows to overlay ground truth points (given in a Nx4 matrix) and the median error of
            %       the fundamental matrix over those points.
            
            %Draw figures
            fig1 = figure; imagesc(self.data(:,:,frameNumber,1)); title('Left image'); hold on;
            posL = get(fig1,'Position');
            set(fig1, 'Position', posL - [posL(1) 0  0 0]/1.2);
            ax1 = gca;
            
            fig2 = figure; imagesc(self.data(:,:,frameNumber,2)); title('Right image'); hold on;
            posR = get(fig2,'Position');
            set(fig2, 'Position', posR + [posR(1) 0  0 0]/1.2);
            ax2 = gca;
            xLim = get(ax2,'xlim')';
            
            %Ground Truth Points
            if exist('groundTruth', 'var')
                GroundTruthPoints = zeros(size(groundTruth,1),2);
                
                for iterPoints = 1:size(groundTruth,1)
                    
                    GroundTruthPoints(iterPoints, 1) = plot(groundTruth(iterPoints,1), groundTruth(iterPoints,2),...
                        [self.colors(mod(iterPoints,length(self.colors))+1)...
                        self.shapes(mod( ceil(iterPoints/length(self.colors)) ,length(self.shapes))+1)],...
                        'Parent', ax1);
                    GroundTruthPoints(iterPoints, 2) = plot(groundTruth(iterPoints,3), groundTruth(iterPoints,4),...
                        [self.colors(mod(iterPoints,length(self.colors))+1)...
                        self.shapes(mod( ceil(iterPoints/length(self.colors)) ,length(self.shapes))+1)],...
                        'Parent', ax2);
                    
                end
                
                P1 = [groundTruth(:,1:2) ones(size(groundTruth,1),1)];
                P2 = [groundTruth(:,3:4) ones(size(groundTruth,1),1)];
                distances12 = CalcGeometricDistances(self.fundMat, P1, P2);
                %distances21 = CalcGeometricDistances(self.fundMat, P2, P1);
                
                medDist12 = median(distances12);
                %medDist21 = median(distances21);
                
                %medDist = medDist12 + medDist21;
                
                % % % %                 FP = self.fundMat*[groundTruth(:,1:2)'; ones(1, size(groundTruth,1))];
                % % % %                 weights = sqrt(sum(FP(1:2,:).^2,1));
                % % % %                 distances = abs(diag(([groundTruth(:,3:4), ones(size(groundTruth,1), 1)]*FP))./weights');
                % % % %                 medDist = median(distances);
                % % % %
                % % % %                 FPt = self.fundMat'*[groundTruth(:,3:4)'; ones(1, size(groundTruth,1))];
                % % % %                 weightst = sqrt(sum(FPt(1:4,:).^2,1));
                % % % %                 distances2 = abs(diag(([groundTruth(:,1:2), ones(size(groundTruth,1), 1)]*FPt))./weightst');
                % % % %                 medDist2 = median(distances2);
                % % % %
                % % % %                 medDist = medDist2 + medDist;
                
                outlierPoints = find(distances12>medDist12);
                numOutliers = length(outlierPoints);
                %                 outlierRectanglesR = NaN(numOutliers, 1);
                %                 outlierRectanglesL = outlierRectanglesR;
                
                %                 for iterOutliers = 1:numOutliers
                %                     outlierRectanglesL(iterOutliers) = RotatedRectangle(ax1, groundTruth(outlierPoints(iterOutliers),1:2)', 30, 30, 45, 'r', ':');
                %                     outlierRectanglesR(iterOutliers) = RotatedRectangle(ax2, groundTruth(outlierPoints(iterOutliers),3:4)', 30, 30, 45, 'r', ':');
                %                 end
                
                xlabel({['Mean FundMat Measure: ' num2str(mean(distances12))];['Median FundMat Measure: ' num2str(medDist12)]; [num2str(numOutliers) ' outliers out of ' num2str(size(groundTruth,1)) '.']}, 'Parent', ax1);
                
            end
            
            while (1)
                axes(ax1); %#ok
                [xL, yL] = ginput(1);
                
                if isempty(xL)
                    break;
                else
                    
                    if exist('pointL','var')
                        delete(pointL);
                    end
                    
                    pointL = plot(xL,yL,'rp', 'Parent', ax1);
                    
                    epipolarLine = self.fundMat*[xL; yL; 1];
                    
                    yLim = -(xLim*epipolarLine(1)+epipolarLine(3))/epipolarLine(2);
                    
                    if exist('lineR','var')
                        delete(lineR);
                        clear lineR;
                    end
                    
                    if exist('pointR','var')
                        delete(pointR);
                        clear pointR;
                    end
                    
                    lineR = plot(xLim, yLim,'r', 'Parent', ax2);
                    
                    while (1)
                        axes(ax2); %#ok
                        
                        [xR, yR] = ginput(1);
                        
                        if isempty(xR)
                            break;
                        else
                            xR = round(xR); yR = round(yR);
                        end
                        
                        if exist('pointR','var')
                            delete(pointR);
                            clear pointR;
                        end
                        pointR = plot(xR,yR,'mp', 'Parent', ax2);
                        
                        Q1 = [xLim(1); yLim(1); 1];
                        Q2 = [xLim(2); yLim(2); 1];
                        P = [xR; yR; 1];
                        
                        dist2Line = norm(cross((Q2-Q1),(P-Q1)))/norm(Q2-Q1);
                        
                        xlabel(['Distance to Epipolar Line: ' num2str(dist2Line) ' pixels']);
                        
                    end
                    
                end
            end
        end
        
        function distances = CalculateEpipolarDistances(self, groundTruth)  %CURRENTLY FOR PAIRS ONLY!!
            % medianDistance = SSTV.CalculateEpipolarDistances(groundTruth)
            % Similar to EpipolarDistances, without the display. Instead of
            % displaying the results graphically, only the errors of
            % the fundamental matrix over the GroundTruth points is calculated and returned.
            
                P1 = [groundTruth(:,1:2) ones(size(groundTruth,1),1)];
                P2 = [groundTruth(:,3:4) ones(size(groundTruth,1),1)];
                distances12 = CalcGeometricDistances(self.fundMat, P1, P2);
                %distances21 = CalcGeometricDistances(self.fundMat, P2, P1);
                
                medDist12 = median(distances12);
                %medDist21 = median(distances21);
                
                %medDist = medDist12 + medDist21;
                
                % % % %                 FP = self.fundMat*[groundTruth(:,1:2)'; ones(1, size(groundTruth,1))];
                % % % %                 weights = sqrt(sum(FP(1:2,:).^2,1));
                % % % %                 distances = abs(diag(([groundTruth(:,3:4), ones(size(groundTruth,1), 1)]*FP))./weights');
                % % % %                 medDist = median(distances);
                % % % %
                % % % %                 FPt = self.fundMat'*[groundTruth(:,3:4)'; ones(1, size(groundTruth,1))];
                % % % %                 weightst = sqrt(sum(FPt(1:4,:).^2,1));
                % % % %                 distances2 = abs(diag(([groundTruth(:,1:2), ones(size(groundTruth,1), 1)]*FPt))./weightst');
                % % % %                 medDist2 = median(distances2);
                % % % %
                % % % %                 medDist = medDist2 + medDist;
                
                outlierPoints = find(distances12>medDist12);
                numOutliers = length(outlierPoints);
                %                 outlierRectanglesR = NaN(numOutliers, 1);
                %                 outlierRectanglesL = outlierRectanglesR;
                
                %                 for iterOutliers = 1:numOutliers
                %                     outlierRectanglesL(iterOutliers) = RotatedRectangle(ax1, groundTruth(outlierPoints(iterOutliers),1:2)', 30, 30, 45, 'r', ':');
                %                     outlierRectanglesR(iterOutliers) = RotatedRectangle(ax2, groundTruth(outlierPoints(iterOutliers),3:4)', 30, 30, 45, 'r', ':');
                %                 end
                
                
                distances = distances12';
                
                            
        end
        
        function RelativeCalibration(self, points1, points2, f0, Cx0, Cy0, varargin)
            %   SSTV.RelativeCalibration(method, points1, points2, f0, Cx0, Cy0) - Performs relative calibration according to the
            %       PTV Point Matching algorithm, setting the K, R and t parameters of the SSTV.
            %       <method> is 'full' for full grid search, 'sparse2full' for sparse to full search, or 'separate' for separate grid search.
            %       <points1> and <points2> are two lists of points (each a Nx2 array of pixel coordinates) to be used for optimization,
            %       1 being for the left camera, and 2 for the right. **CURRENTLY ONLY FOR 2 CAMERAS**
            
            if (nargin == 6) %method undefined
                method = 'separate'; %default method - separate.
            
            elseif (nargin == 7) %method defined
                method = varargin{1}; 
                
            else
                disp('Wrong number of input variables');
                keyboard;
            end
            
            if ~( strcmpi(method ,'full') || strcmpi(method ,'sparsefull') || strcmpi(method ,'separate'))
                disp('Wrong method at input');
                keyboard;
            end
            
            %**NO CHECK FOR SAME LENGTH OF POINTS ETC. ASSUME CORRECT INPUT
            
            sigmaDmin = inf; %minimization parameter
            
            Fcurr = self.fundMat;
            
            if (strcmpi(method ,'full'))
                %% Full Grid Search
                
                for iterf = 1:100
                    f = f0/2 + f0/100*iterf;
                    for iterCx = 1:30
                        Cx = Cx0-15 + iterCx;
                        
                        for iterCy = 1:30
                            Cy = Cy0-15 + iterCy;
                            
                            Kcurr = self.CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K
                            
                            % tic;
                            [Rcurr, tcurr, Fcurr, sigmaD] = calcRtFromFK(Fcurr, Kcurr, points1, points2);
                            
                            if sigmaD < sigmaDmin %improvement in sigmaD
                                disp (['[f ' num2str(f) ', Cx' num2str(Cx) ', Cy' num2str(Cy) ']: sigmaD = ' num2str(sigmaD)]);
                                
                                sigmaDmin = sigmaD;
                                RBest = Rcurr;
                                tBest = tcurr;
                                
                                FBest = Fcurr;
                                
                                fBest = f;
                                CxBest = Cx;
                                CyBest = Cy;
                                
                                KBest = CalcKFromParameters(fBest, CxBest, CyBest);
                            end
                            % toc
                        end
                        
                    end
                    
                    %keepalive
                    disp('');
                    disp(['Iteration ' num2str(iterf) 'in f done']);
                end
                
            elseif (strcmpi(method ,'sparse2full'))
                %% Sparse2Full Grid Search
                %tic;
                for iterf = 0:10
                    f = f0/2 + f0/10*iterf;
                    for iterCx = 0:2
                        Cx = Cx0-15 + iterCx*15;
                        
                        for iterCy = 0:2
                            Cy = Cy0-15 + iterCy*15;
                            
                            Kcurr = CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K
                            
                            % tic;
                            [Rcurr, tcurr, Fcurr, sigmaD] = calcRtFromFK(Fcurr, Kcurr, points1, points2);
                            
                            if sigmaD < sigmaDmin %improvement in sigmaD
                                disp (['[f ' num2str(f) ', Cx' num2str(Cx) ', Cy' num2str(Cy) ']: sigmaD = ' num2str(sigmaDS2F)]);
                                
                                sigmaDmin = sigmaD;
                                RBest = Rcurr;
                                tBest = tcurr;
                                
                                FBest = Fcurr;
                                
                                fBest = f;
                                CxBest = Cx;
                                CyBest = Cy;
                                
                                KBest = CalcKFromParameters(fBest, CxBest, CyBest);
                            end
                            % toc
                        end
                        
                    end
                    
                    disp('');
                    disp(['Iteration ' num2str(iterf) 'in f done']);
                    
                end
                toc;
                disp('');
                disp('Starting Refinement');
                
                %tic;
                for iterf = 0:10
                    f = fBest-f0/100*5 + f0/100*iterf;
                    for iterCx = 0:3
                        Cx = CxBest-5 + iterCx*5;
                        
                        for iterCy = 0:3
                            Cy = CyBest-5 + iterCy*5;
                            
                            Kcurr = CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K
                            
                            % tic;
                            [Rcurr, tcurr, Fcurr, sigmaD] = calcRtFromFK(Fcurr, Kcurr, points1, points2);
                            
                            if sigmaD < sigmaDmin %improvement in sigmaD
                                disp (['[f ' num2str(f) ', Cx' num2str(Cx) ', Cy' num2str(Cy) ']: sigmaD = ' num2str(sigmaDS2F)]);
                                
                                sigmaDmin = sigmaD;
                                RBest = Rcurr;
                                tBest = tcurr;
                                
                                FBest = Fcurr;
                                
                                fBest = f;
                                CxBest = Cx;
                                CyBest = Cy;
                                
                                KBest = CalcKFromParameters(fBest, CxBest, CyBest);
                            end
                            % toc
                        end
                        
                    end
                    
                    disp('');
                    disp(['Iteration ' num2str(iterf) 'in f (refined) done']);
                    
                end
            else %Separate Grids
                %tic;
                Cx = Cx0;
                Cy = Cy0;
                for iterf = 1:100
                    f = f0*3/4 + f0/100*iterf;
                    
                    Kcurr = CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K
                    
                    [Rcurr, tcurr, Fcurr, sigmaD] = calcRtFromFK(Fcurr, Kcurr, points1, points2);
                    
                    if sigmaD < sigmaDmin %improvement in sigmaD
                        disp (['[f ' num2str(f) ', Cx' num2str(Cx) ', Cy' num2str(Cy) ']: sigmaD = ' num2str(sigmaD)]);
                        
                        sigmaDmin = sigmaD;
                        RBest = Rcurr;
                        tBest = tcurr;
                        
                        FBest = Fcurr;
                        
                        fBest = f;
                        CxBest = Cx;
                        CyBest = Cy;
                        
                        KBest = CalcKFromParameters(fBest, CxBest, CyBest);
                    end
                    
                end
                
                disp('');
                disp('Iterations in f done');
                
                f = fBest;
                for iterCx = 1:30
                    Cx = Cx0-15 + iterCx;
                    
                    Kcurr = CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K
                    
                    [Rcurr, tcurr, Fcurr, sigmaD] = calcRtFromFK(Fcurr, Kcurr, points1, points2);
                    
                    if sigmaD < sigmaDmin %improvement in sigmaD
                        disp (['[f ' num2str(f) ', Cx' num2str(Cx) ', Cy' num2str(Cy) ']: sigmaD = ' num2str(sigmaD)]);
                        
                        sigmaDmin = sigmaD;
                        RBest = Rcurr;
                        tBest = tcurr;
                        
                        FBest = Fcurr;
                        
                        fBest = f;
                        CxBest = Cx;
                        CyBest = Cy;
                        
                        KBest = CalcKFromParameters(fBest, CxBest, CyBest);
                    end
                    
                    
                end
                disp('');
                disp('Iterations in Cx done');
                
                f = fBest;
                Cx = CxBest;
                for iterCy = 1:30
                    Cy = Cy0-15 + iterCy;
                    
                    
                    Kcurr = CalcKFromParameters(f, Cx, Cy); %K1 = K2 = K
                    
                    [Rcurr, tcurr, Fcurr, sigmaD] = calcRtFromFK(Fcurr, Kcurr, points1, points2);
                    
                    if sigmaD < sigmaDmin %improvement in sigmaD
                        disp (['[f ' num2str(f) ', Cx' num2str(Cx) ', Cy' num2str(Cy) ']: sigmaD = ' num2str(sigmaD)]);
                        
                        sigmaDmin = sigmaD;
                        RBest = Rcurr;
                        tBest = tcurr;
                        
                        FBest = Fcurr;
                        
                        fBest = f;
                        CxBest = Cx;
                        CyBest = Cy;
                        
                        KBest = CalcKFromParameters(fBest, CxBest, CyBest);
                    end
                    % toc
                end
                disp('');
                disp('Iterations in Cy done');
                
                disp('');
                disp('all iterations complete.');
                toc;
            end
            
            self.fundMat = FBest;
            self.K = KBest;
            self.R = RBest;
            self.t = tBest;
                
        end
        
        
        %OVERLOADED METHODS
        function varagout = subsref(self,S)
            switch S(1).type
                case '()'
                    if (length(S) < 2)
                        varagout = builtin('subsref', self.data,S);
                        return;
                    else
                        varagout = builtin('subsref',self,S);
                    end
                    
                case '.'
                                        if nargout
                                            varagout = builtin('subsref',self,S);
                                        else
                                            builtin('subsref',self,S);
                                        end
%                     varagout = builtin('subsref',self,S);
                case '{}'
                    disp('{} is Not a supported subscripted reference for SSTV')
                    keyboard;
            end
        end
        
        function sizeOut = size(self,varargin)
            if nargin > 1
                sizeOut = size(self.data,varargin{:});
            else
                sizeOut = size(self.data);
            end
        end
        
    end % methods
end % class