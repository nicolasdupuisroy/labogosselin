function out = imEnergy(im,area)

% imEnergy    This function computes the energy of the signal (image in
% this case). The definition is taken from Jason Gold PhD dissertation
% (Signal and Noise in Perceptual Learning).  
% 
%     imEnergy(im,area) Computes the root mean square contrast of the image
%     squared and multiplied by the image area in degree^2 of visual angle.
% 
% See also RMScontrast, noisePowerSpectralDensity
% 
% Nicolas Dupuis-Roy, 2012-10-10

rmsc = RMScontrast(im);
out = rmsc^2 * area;

