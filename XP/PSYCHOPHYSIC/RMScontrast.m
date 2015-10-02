function varargout = RMScontrast(in)

% RMScontrast     This function computes two definition of the root mean square contrast. 
% 
%     First one is described in Jason Gold PhD dissertation (Signal and
%     Noise in Perceptual Learning). The second one is described in
%     Wikipedia at
%     http://en.wikipedia.org/wiki/Contrast_(vision)#RMS_contrast
% 
% See also noisePowerSpectralDensity, imEnergy
% 
% Nicolas Dupuis-Roy, 2012-10-10

in = double(in);
if max(in)>1
    in = in/255;
end

L = mean(in(:));
c = (in(:)-L)/L;
varargout{1} = sqrt(sum(c.^2)/numel(c));
 
if nargout==2
    varargout{1} = sqrt(sum((in(:)-mean(in(:))).^2)/numel(in));
end

