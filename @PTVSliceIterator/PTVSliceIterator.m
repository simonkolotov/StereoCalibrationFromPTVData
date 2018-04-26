classdef PTVSliceIterator < handle
    %PTVSliceIterator is a class for the iterator of slices in the PTV Point Matching Algorithm.
    %   The iterator defines the frames to be used in slicing the data, as well as the frames to be used for retrieving the points from the matched tracks. It
    %   also fascilitates the process of slicing detection, etc.
    %
    %FIELDS:
    %   EndOfSequence: boolean (1 or 0) indicating end of sequence (no more slices)
    %   
    %   firstFrameOfSlice, lastFrameOfSlice, FramesPerSlice, lastFrameOfSequence: as it sounds   ***Change here for a different implementation***
    %
    %CONSTRUCTORS:
    %   PTVSliceIteratorInstance = PTVSliceIterator() - Create empty object of PTVSliceIterator Class
    %
    %   PTVSliceIteratorInstance = PTVSliceIterator(PTVSIInstanceOld) - Initialize from another PTVSI Instance <PTVSIInstanceOld>
    %
    %   PTVSliceIteratorInstance = PTVSliceIterator(superSpaceTimeVolumeIn) - Initialize by taking the number of frames in the sequence represented by the SSTV,
    %       and using the default number of frames per slice (15).
    %       
    %   PTVSliceIteratorInstance = PTVSliceIterator(superSpaceTimeVolumeIn, framesPerSlcie) - Initialize by taking the number of frames in the sequence represented by the SSTV,
    %       and using the given number of frames per slice.
    %
    %METHODS:
    %   PTVSI.NextSlice() - move to the next slice by adding framesPerSlice to the first and the last frame of the slice.
    %       if lastFrameOfSlice moves beyond the sequence limits, the EndOfSequence flag is lifted.
    %   
    %Simon Kolotov, Ver 2.0, Spring 2013
    
    properties (SetAccess = public, GetAccess = public)
        
        EndOfSequence = 0;
        
        framesPerSlice = 15;
        
        firstFrameOfSlice = 1;
        
        lastFrameOfSlice
        
        lastFrameOfSequence
        
    end %public Properties
    
    methods
        %CONSTRUCTORS
        %   PTVSliceIteratorInstance = PTVSliceIterator() - Create empty object of PTVSliceIterator Class
        %
        %   PTVSliceIteratorInstance = PTVSliceIterator(PTVSIInstanceOld) - Initialize from another PTVSI Instance <PTVSIInstanceOld>
        %
        %   PTVSliceIteratorInstance = PTVSliceIterator(superSpaceTimeVolumeIn) - Initialize by taking the number of frames in the sequence represented by the SSTV,
        %       and using the default number of frames per slice (15).
        %
        %   PTVSliceIteratorInstance = PTVSliceIterator(superSpaceTimeVolumeIn, framesPerSlcie) - Initialize by taking the number of frames in the sequence represented by the SSTV,
        %       and using the given number of frames per slice.
        function PTVSliceIteratorInstance = PTVSliceIterator(varargin)
            
            PTVSliceIteratorInstance.lastFrameOfSlice = PTVSliceIteratorInstance.framesPerSlice;
            
            %   ***
            if (nargin < 1)
                PTVSliceIteratorInstance.lastFrameOfSequence = 0;
            
            elseif (nargin == 1)
                
                if ( strcmpi(class(varargin{1}) ,'PTVSliceIterator') ) %Initialize from another PTVSI
                    PTVSliceIteratorInstance.EndOfSequence = varargin{1}.EndOfSequence;
                    PTVSliceIteratorInstance.framesPerSlice = varargin{1}.framesPerSlice;
                    PTVSliceIteratorInstance.lastFrameOfSlice = varargin{1}.lastFrameOfSlice;
                    PTVSliceIteratorInstance.EndOfSequence = varargin{1}.EndOfSequence;
                    PTVSliceIteratorInstance.lastFrameOfSequence = varargin{1}.lastFrameOfSequence;
                     
                elseif ( strcmpi(class(varargin{1}) ,'SuperSpaceTimeVolume') ) %Initialize from SSTV
                    PTVSliceIteratorInstance.lastFrameOfSequence = size(varargin{1}, 3);
                    
                else
                    disp('Wrong Class at input of PTVSI Constructor');
                    keyboard;
                end
                
            elseif (nargin == 2) %Initialize from SSTV and numFrames
                if ( strcmpi(class(varargin{1}) ,'SuperSpaceTimeVolume') ) %Initialize from SSTV
                    PTVSliceIteratorInstance.lastFrameOfSequence = size(varargin{1}, 3);
                    PTVSliceIteratorInstance.framesPerSlice = varargin{2};
                end
                
            else
                disp('Wrong Number of parameters at input of PTVSI Constructor');
                keyboard;
                
            end %nargin check
        end %constructor
        
        %METHODS
        function self = NextSlice(self)
            %   PTVSI.NextSlice() - move to the next slice by adding framesPerSlice to the first and the last frame of the slice.
            %       if lastFrameOfSlice moves beyond the sequence limits, the EndOfSequence flag is lifted.
            
            self.firstFrameOfSlice = self.firstFrameOfSlice + self.framesPerSlice;
            
            self.lastFrameOfSlice = self.lastFrameOfSlice + self.framesPerSlice;
            
            if ( self.lastFrameOfSlice > self.lastFrameOfSequence)
                
                self.EndOfSequence = 1;
                
            end
            
        end
               
    end %Public Methods
    
end % class