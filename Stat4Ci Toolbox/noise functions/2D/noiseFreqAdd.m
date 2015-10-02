function varargout = noiseFreqAdd(varargin);% output = noiseFreqAdd([xSize,ySize], noiseFunction, prop, 'add', snr, profile)% ------------------------------------------------------------------------% -- 0.0 global and local variables initializationglobal glopar;if ((length(glopar) < 4) | (length(glopar) > 7) )	s1 = sprintf('glopar{1}[x, y] = size along the x and y axis\nglopar{2} the pdf to be chosen (e.g. gaussian or uniform)\n');	s2 = sprintf('glopar{3} proportion of active pixel  (1 = all pixels active)glopar{4} = signal-to-noise ratio \n');	s3 = sprintf('glopar{5} the mixing funtion (e.g. add, mult or speckle)\nglopar{6} = a vector (0 = no filtering), the filter profile \n');	s4 = sprintf('Must be defined in the CID.\nglopar{7} = information about the stimulus: bits, min and max\n');	error([s1,s2,s3,s4]);end% global variablexSize = glopar{1}(1);ySize = glopar{1}(2);noiseFunction = glopar{2}prop = glopar{3}noiseControl = glopar{4};if length(glopar) < 5,	mixFunction = 'add';else,					mixFunction = glopar{5}end;if length(glopar) < 6,		filteredNoise = 0;else	if glopar{6} == 0		filteredNoise = 0;	else		profile = glopar{6};		filteredNoise = 1	end;end;% inputoneTrial = varargin{1};if xSize ~= ySize,	error('the images must be squared');end;	% ------------------------------------------------------------------------% -- 1.0 generate the noiseswitch noiseFunctioncase 'gaussian'	randn('state',  oneTrial(1));	rand('state',  oneTrial(1));	aWhiteNoise = randn(xSize, ySize) .* round(rand(xSize, ySize)-(.5-prop));case 'uniform'	rand('state',  oneTrial(1));	aWhiteNoise = rand(xSize, ySize) .* round(rand(xSize, ySize)-(.5-prop));otherwise	error('the noise generating function is not specified')end;muTheNoise = mean(aWhiteNoise(:));% ------------------------------------------------------------------------% -- 2.0 filter the noise spectrumif filteredNoise	[x, y] = meshgrid(-xSize/2:xSize/2-1, -ySize/2:ySize/2-1);	rayon = sqrt(x .^2 + y .^2);	rayon = max(min(rayon, round(xSize / 2)), 1);	noise_f = fft2(aWhiteNoise);	noise_f_filt = noise_f .* ifftshift(profile(round(rayon)));	theNoise = real(ifft2(noise_f_filt));	muTheNoise = mean(theNoise(:));else	theNoise = aWhiteNoise;end;if nargin == 1	varargout{1} = theNoise;	return;end;% ------------------------------------------------------------------------% -- 3.0 mix the noise with an imageanImage = double(varargin{2});anImage = (anImage - min(anImage(:))) / (max(anImage(:)) - min(anImage(:)));switch mixFunctioncase 'add'	nrj_Im = sum(sum( (anImage-.5).^2 ));	nrj_Noi = sum(sum( (theNoise-muTheNoise).^2 ));	S_N = 10*log10(nrj_Im/nrj_Noi)	theStimulus = anImage + theNoise / noiseControl * S_N;case 'mult'	% semble inutile !!!!!!!!	theStimulus =  noiseControl * anImage .* theNoise;case 'speckle'	theStimulus = anImage + sqrt(12*noiseControl) * anImage .* theNoise; % varNoise = 0.04end;% ------------------------------------------------------------------------% -- 4.0 Output and infornation on the stimulusmmin = min(theStimulus(:));mmax = max(theStimulus(:)); % we suppose that (1) images are scaled beetween 0 and 1 and (2) we used 256 gray level, i.e. max bits = 8bits = log2(floor((mmax-mmin)*256));glopar{7}(1) = bits;glopar{7}(2) = mmin;glopar{7}(3) = mmax;varargout{1} = theStimulus;