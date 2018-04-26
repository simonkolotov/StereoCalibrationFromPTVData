classdef SuperSlice < handle
    %superSlice is a class defining a "slice" in time of super space-time volume for the PTV Point Matching Algorithm
    %
    %DATA ACCESS:
    %   SS(X, Y, CamView), where:
    %                   (X,Y) are the coordiantes in a single Slice Image,
    %                   CamView is the Index of the camera in the Camera Array
    %
    %CONSTRUCTORS:
    %   superSliceInstance = SuperSlice() - Create empty object of superSlice Class
    %
    %   superSliceInstance = SuperSlice(sSInstanceOld) - Initialize from another SS Instance <sSInstanceOld>
    %      
    %   superSliceInstance = SuperSlice(sSTVInstance, frameGroup) - Initialize from a slice of an <SSTVInstance> defined by <frameGroup>,
    %       which is a vector of frame indices (e.g. 1:15).
    %
    %METHODS:
    %   SS.GetView(cameraNum) - retrieve the slice view from camera <cameraNum>
    %
    %Simon Kolotov, Ver 2.0, Spring 2013
    
    properties (SetAccess = public, GetAccess = public)
        data
    end %Public Properties
    
    methods
        %Constructors
        function superSliceInstance = SuperSlice(varargin)
            %   superSliceInstance = superSlice() - Create empty object of superSlice Class
            %
            %   superSliceInstance = superSlice(sSInstanceOld) - Initialize from another SS Instance <sSInstanceOld>
            %
            %   superSliceInstance = superSlice(sSTVInstance, firstFrame, lastFrame) - Initialize from a slice of an <SSTVInstance> defined by the limits <firstFrame>
            %       and <lastFrame>
            
            if (nargin==1) %Initialize from another SS
                
                if (strcmpi(class(varargin{1}) ,'SuperSlice')) 
                    superSliceInstance.data = varargin{1}.data;
                else
                   disp('Wrong Class at input of SS Constructor');
                   keyboard;
                end
                
            elseif (nargin==2) %Initialize from SSTV
                
                if ~(strcmpi(class(varargin{1}) ,'SuperSpaceTimeVolume')) 
                    disp('Wrong Class at input of SS Constructor');
                    keyboard;
                end
                
                frameGroup = varargin{2};
                
                superSliceInstance.data = squeeze(max(varargin{1}.data(:,:,frameGroup, :),[],3));
                
            elseif (nargin > 2)
                disp('Wrong number of parameters at input of STL Constructor');
                keyboard;
            
            end  %nargin check
        
        end %Constructor
        
        function view = GetView(self, cameraNum)
            %SSTV.GetView(cameraNum)
            %
            %retrieve the slcie view from camera <cameraNum>
            view = self.data(:,:,cameraNum);
        end
        
        %OVERLOADED METHODS
        function data = subsref(self,S)
            switch S(1).type
                case '()'
                    if (length(S) < 2)
                        data = builtin('subsref', self.data,S);
                        return;
                    else
                        data = builtin('subsref',self,S);
                    end
                    
                case '.'
                    data = builtin('subsref',self,S);
                    
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