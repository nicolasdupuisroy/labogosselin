function out = overlay_gradual(varargin)

%overlay_gradual    This function is mainly for exhibition purposes. In this
%     case it multiplies (or adds) a classification image with (to) a template
%     image. 
% 
%     out = overlay_gradual(im,zscoresMap,nLevels,maxContrast) shows
%     nLevels (zscoresMap*maxContrast X image). The nLevels go from
%     -maxContrast to +maxContrast. The type of overlay is multiplicative
%     by default. You can also add a fifth argument specifying the type,
%     either 'multiplicative' or 'additive'
%         
% 
% Nicolas Dupuis-Roy
% 2012-10-01
% 
% See also overlay_cluster, overlay_pixel and find_sig_cluster

if nargin<4
    error('This function needs at least 4 input arguments')
end

% Assign input argument to variable
IM = varargin{1};
zscores = varargin{2};
nLevels = varargin{3};
maxContrast = varargin{4};
if nargin<5
    type = 'multiplicative';
else
    type = varargin{5};
end

zscores = stretch(sameSize(IM,zscores));
zweight = linspace(-maxContrast,maxContrast,nLevels);
sign = (zweight>0) *2-1;
out = zeros(size(IM,1),size(IM,2),size(IM,3), nLevels);
% IM = stretch(IM);
switch type
    case 'additive'
        for ii = 1:nLevels
            tmp = zscores;
            tmp = stretch(sign(ii)*tmp);
            out(:,:,:,ii) = IM*(1-abs(zweight(ii))) + tmp*abs(zweight(ii));
        end
    case 'multiplicative'
        for ii = 1:nLevels
            tmp = zscores;
            tmp = stretch(sign(ii)*tmp);
            out(:,:,:,ii) = IM .* (tmp*abs(zweight(ii))+(1-abs(zweight(ii))));
        end
end        
if nargin==5
    figure,
    for ii = 1:nLevels 
        subplot(1,nLevels,ii), imshow(out(:,:,:,ii))
    end
end
