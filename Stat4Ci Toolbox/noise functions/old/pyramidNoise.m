function varargout=pyramidNoise(varargin)global gloparif (length(glopar)~=4)	error(sprintf('glopar{1} = size along the X axis\nglopar{2} = size along the Y axis\nglopar{3} = std of one bubble in minimum number of cycles\nglopar{4} = a vector of 5 numbers of bubbles\nMust be defined in the CID.'));endnBands=6;xSize = glopar{1};ySize = glopar{2};oneTrial = varargin{1};rand('state', oneTrial(1));theNoise=zeros(ySize,xSize,nBands-1);for ii=1:(nBands-1),	nn=0;	tempRand = rand(ySize,xSize);	for jj=1:glopar{4}(ii),		[scrap xMax]=max(max(tempRand));		[scrap yMax]=max(tempRand(:,xMax));  		theNoise(yMax,xMax,ii)=theNoise(yMax,xMax,ii)+1;		tempRand(yMax,xMax)=0;	endendvarargout{1} = theNoise;if (nargin == 2)	anImage = varargin{2};	nBands=6;	[ySize, xSize] = size(anImage);	nPeriod=glopar{3}; nZero=2.18;	stdev = nPeriod * 2^(5);	maxHalfSize5 = round(stdev * nZero);	gauss5 = zeros(2*maxHalfSize5,2*maxHalfSize5);	[y,x] = meshgrid(-maxHalfSize5:maxHalfSize5,-maxHalfSize5:maxHalfSize5);	gauss5 = exp(-(x.^2/stdev^2)-(y.^2/stdev^2));	gauss5 = gauss5/max(gauss5(:));	clear x, y;	method = 'nearest';	stdev = nPeriod * 2^(4);	maxHalfSize4 = round(stdev * nZero);	gauss4 = double(imresize(gauss5,[(2*maxHalfSize4+1) (2*maxHalfSize4+1)], method));	stdev = nPeriod * 2^(3);	maxHalfSize3 = round(stdev * nZero);	gauss3 = double(imresize(gauss5,[(2*maxHalfSize3+1) (2*maxHalfSize3+1)], method));	stdev = nPeriod * 2^(2);	maxHalfSize2 = round(stdev * nZero);	gauss2 = double(imresize(gauss5,[(2*maxHalfSize2+1) (2*maxHalfSize2+1)], method));	stdev = nPeriod * 2^(1);	maxHalfSize1 = round(stdev * nZero);	gauss1 = double(imresize(gauss5,[(2*maxHalfSize1+1) (2*maxHalfSize1+1)], method));	fGauss=mkGaussian(11);	fGauss1=fGauss(:,1);	fGauss1=sqrt(2)*fGauss1/sum(fGauss1);	clear fGauss;	winPlane = double(zeros(ySize,xSize));	[pyr,pind] = buildLpyr(double(anImage),nBands,fGauss1);	stimulus = zeros(ySize,xSize);	for ii = 1:(nBands-1),			nameGauss = sprintf('gauss%d',ii);		nameMax = sprintf('maxHalfSize%d',ii);		temp = eval(nameMax);		tempPlane = zeros(ySize+temp-1,xSize+temp-1);		tempPlane = real(ifft2(fft2(eval(nameGauss),ySize+temp-1,xSize+temp-1).*fft2(theNoise(:,:,ii), ySize+temp-1, xSize+temp-1)));		winPlane = min(tempPlane(temp:ySize+temp-1,temp:xSize+temp-1), 1);		stimulus = stimulus + double(reconLpyr(pyr,pind,[ii],fGauss1).*winPlane);	end	varargout{2} = uint8(stimulus + double(reconLpyr(pyr,pind,[nBands],fGauss1)));end