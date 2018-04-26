function RenameData
% Renames image data collected at Alex Liberzon's lab from lab standard
% names to Left_1/Right_1 and Left_2/Right_2 pairs.

close all;
disp(' '); disp(' ');
home;

%Function Settings
askUserForPath = 0;

%Default Data Directory
dataPath = 'D:\Personal\Senya\Dropbox\Thesis\Data\Dataset 2_03\';

%Default Input Directory Addendum
% pathIn = '\Raw\Img';
pathIn = '';
%Default Output Directory Addendum
pathOut = '\Images';

%Resultant image format
outImgFormat = '.tif';

%Number of cameras
numCams = 4;



%% Select Data folder and check its integrity
%Select root directory
if askUserForPath
    %Set Path to Data
    dataPathUser = uigetdir(dataPath,...
        'Select the Data Directory, Cancel to Quit');
    if dataPathUser
        dataPath = dataPathUser;
    else
        home;
        disp('A valid path is required to load data!');
        return
    end
end

%Locate Raw Image Directory
rawImageDir = [dataPath pathIn];
if ~isdir(rawImageDir)
    %No Raw Images directory exists at default place - ask user.
    disp('Raw Images Directory not found... Where could it be?..');
    dataPathUser = uigetdir(dataPath,...
        'Select the Raw Image Directory, Cancel to Quit');
    if dataPathUser
        rawImageDir = dataPathUser;
        disp('Is that it?...');
    else
        disp('A valid path is required to the raw image data!');
        return;
    end
else
    disp('Raw Image Directory Located.');
end

%Locate Raw Image Files
disp('Enumerating Raw Files...');
rawFileNames = dir([rawImageDir '\cam*']);

%Check for files in Raw Images directory
if isempty(rawFileNames)
    disp('No Raw Files Found!');
    return
    %Check 4 cameras,_target file for each image (BASIC CHECK - IMPROVE IF NECESSARY)
elseif mod(length(rawFileNames),2*numCams)
    disp('Wrong file quantity - check integrity of raw data!');
    return
else    
    disp('OK.')
    
    rawFileNames = rawFileNames(1:2:end); %remove _target files from list
    
    numFiles = length(rawFileNames);
    numFrames = numFiles/numCams;
    orderMagOfFrames = 10^(length(num2str(numFrames))-1);   %Which order of magnitde is the number of frames (For trailing zeros in names)
end

%Locate Final Image Directory
finalImageDir = [dataPath pathOut];
if ~isdir(finalImageDir)
    %No Final Images directory exists at default place - ask user.
    disp('Output Images Directory not found... Where could it be?..');
    dataPathUser = uigetdir(dataPath,...
        'Select the Raw Image Directory, Cancel to Quit');
    if dataPathUser
        finalImageDir = dataPathUser;
        disp('Is that it?...');
    else
        disp('A valid path is required to save final images!');
        return
    end
else
    disp('Output Image Directory Located.');
end

%Check for existing Image Files
disp('Checking for files in Output Image Directory...');
outputFiles = dir(finalImageDir);
if ~isdir(outputFiles(end).name)
    %there are files in the target directory
    overwriteFiles = abs(menu('Target Directory Occupied!', 'Overwrite', 'No, let me check first!') - 2);
    if ~overwriteFiles
        return;
    end
else
    disp('OK.')
end


%Copy the files
disp('Copying the files...');
for iterRawFile = 1:numFiles
    %Prefix to file
    switch ceil(iterRawFile/numFrames)           %CASE TAILORED FOR 4 CAMS!
        case 1  %cam1
            sidePrefix = '\Left_1 ';
        case 2  %cam2
            sidePrefix = '\Right_1 ';
        case 3  %cam3
            sidePrefix = '\Right_2 ';
        case 4  %cam4
            sidePrefix = '\Left_2 ';
    end
    
    %Number of trailing zeros
    frameNumber = mod(iterRawFile,numFrames);
    if frameNumber == 0
        frameNumber = numFrames;
    end
    
    trailingZeros = ceil(log10(orderMagOfFrames/frameNumber));
    
    numberPrefix = repmat('0',[1 trailingZeros]);
    
    %Copy the file
    copyStatus = copyfile(...
        [rawImageDir '\' rawFileNames(iterRawFile).name],...
        [finalImageDir, sidePrefix, numberPrefix, num2str(frameNumber), outImgFormat]);
        
    if ~copyStatus %Some error
        error(['Problem in file ' num2str(iterRawFile)]);
    end
    
    if ~mod(iterRawFile,100)  %Still Alive and Kicking
        disp(['File ' num2str(iterRawFile) ' Done']);
    end
    
end

disp('All Files Copied.')
end