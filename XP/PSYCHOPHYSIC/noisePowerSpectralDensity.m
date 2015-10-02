function out = noisePowerSpectralDensity(noise,area)

% noisePowerSpectralDensity   This function computes the noise power
% spectral density as described at page 34 of Jason Gold PhD dissertation
% (Signal and Noise in Perceptual Learning).   
% 
%     noisePowerSpectralDensity(noise,areaOfSinglePixel) The area of the
%     image must be in degrees^2 of visual angle.
% 
% See also imEnergy, RMScontrast
% 
% Nicolas Dupuis-Roy, 2012-10-10

out = std(noise(:))^2 * (area/numel(noise));