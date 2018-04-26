function main
%the function configures and runs the PTV Point Matching Algorithm.
%
%Simon Kolotov, Ver 2.0, Spring 2013

% clear classes;
close all;

downsampleRate = 4;

percentageOfBestPoints = .5; %How many of the points to keep for final parametrization. 0.5 - to fit the median

%Initial Calibration Parameters
f0 = 100/.017; Cx0 = 512; Cy0 = 640; %focal length/nPixelsPerMm, half width and half height

% pathIn = 'D:\Simon Kolotov\Dropbox\Thesis\Data\14 Longish Datasets\Dataset 2_03\';
pathIn = 'D:\Personal\Senya\Dropbox\Thesis\Data\14 Longish Datasets\Dataset 2_06\';
names = {'Left_'; 'Right_'};

%% Load data into SSTV

tic;
SSTV = SuperSpaceTimeVolume(pathIn, names);  %ca. 9.5 sec
SSTV_Downsampled = SSTV.Downsample(downsampleRate);       %ca. 2.5 sec

SliceIterator = PTVSliceIterator(SSTV_Downsampled);

pointList = CorrespondingPointList();

tac('Loading data')
%% Iterate through slices
tic;
while ~SliceIterator.EndOfSequence          %ca. .3 sec per slice
    %Create Slice
    slice = SuperSlice(SSTV_Downsampled, ( SliceIterator.firstFrameOfSlice : SliceIterator.lastFrameOfSlice ) );
    
    %Detect Tracks + Match Tracks
    trackList = SuperTrackList(slice);
    matchedTracks = trackList.Match();
    
    %Detect Points
    pointList.AddPointsFromAnotherSlice(matchedTracks, SSTV_Downsampled, [SliceIterator.firstFrameOfSlice SliceIterator.lastFrameOfSlice]);

    SliceIterator.NextSlice();
end %ca. 15 sec for 10 slices
tac('Iterations');
%% Fundamental Matrix
disp('Calculating Fundmat...');
tic;
SSTV_Downsampled.fundMat = pointList.CalcFundMat('IRLS');
tac('DS');
SSTV.fundMat = SSTV_Downsampled.ResampleF(1/downsampleRate);
tac('Full');
disp('Done');

% load([pathIn 'myGroundTruth Dataset 2_02.mat']);

%%Check Fundmat Results
%SSTV.EpipolarDistances(1, pointList.data*downsampleRate);


%% Select the best (Highest weight) points, and calibrate.
%Best 50% points (up to median)
[~, sortedLocations] = sort(pointList.pointWeights, 1, 'descend');

nBestPoints = size(pointList.data,1)*percentageOfBestPoints;
points = pointList.data(sortedLocations(1:nBestPoints),:)*downsampleRate;

points1 = [points(:,1:2)'; ones(1, nBestPoints)];
points2 = [points(:,3:4)'; ones(1, nBestPoints)];

disp('Calculating the calibration...')  ;
tic;
SSTV.RelativeCalibration(points1, points2, f0, Cx0, Cy0, 'separate');
tac('Separate Calibration');
disp('Done.');
end