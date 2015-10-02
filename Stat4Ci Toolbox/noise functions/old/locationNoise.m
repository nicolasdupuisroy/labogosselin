function varargout=locationNoise(varargin)global gloparif (length(glopar)~=4)	error(sprintf('glopar{1} = size along the X axis\nglopar{2} = size along the Y axis\nglopar{3} = std of one bubble\nglopar{4} = number of bubbles\nMust be defined in the CID.'));endoneTrial = varargin{1};rand('state', oneTrial(1));xSize = glopar{1};ySize = glopar{2};theNoise = zeros(ySize, xSize);tempRand = rand(ySize, xSize);for jj = 1:glopar{4},	[scrap xMax]=max(max(tempRand));	[scrap yMax]=max(tempRand(:,xMax));	theNoise(yMax,xMax)=theNoise(yMax,xMax)+1;	tempRand(yMax,xMax)=0;endvarargout{1} = theNoise;if (nargin == 2)	anImage = varargin{2};	stdev = glopar{3};	nZero = 3;	maxHalfSize = round(stdev * nZero);	gauss = zeros(2*maxHalfSize,2*maxHalfSize);	[y,x] = meshgrid(-maxHalfSize:maxHalfSize,-maxHalfSize:maxHalfSize);	gauss = exp(-(x.^2/stdev^2)-(y.^2/stdev^2));	gauss = gauss/max(gauss(:));	clear x, y;	winPlane = zeros(ySize,xSize);	stimulus = zeros(ySize,xSize);	tempPlane = zeros(ySize+maxHalfSize-1,xSize+maxHalfSize-1);	tempPlane = real(ifft2(fft2(gauss,ySize+maxHalfSize-1,xSize+maxHalfSize-1).*fft2(theNoise, ySize+maxHalfSize-1, xSize+maxHalfSize-1)));	winPlane = min(tempPlane(maxHalfSize:ySize+maxHalfSize-1,maxHalfSize:xSize+maxHalfSize-1), 1);	stimulus = uint8(round(winPlane.*(double(anImage)-128))+128);	varargout{2} = stimulus;elseif (nargin >2)	error('Too many input arguments.')end