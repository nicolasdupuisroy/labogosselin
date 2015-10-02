function varargout = sameSize(varargin)

%sameSize    Resize all matrix so that they have the same x and y size.
%   [out1 out2] = sameSize(in,in2) Change the size of in2 so that it
%   matches the x- and y-size of in. This functions takes as many arguments
%   as there are images.
% 
% 
% Nicolas Dupuis-Roy
% 2012-10-01
%
%   See also overlay, overlay_pixel and overlay_cluster

if (nargin-1)~=nargout
    help sameSize
    error('wrong number of input/output arguments')
end

im = varargin{1};
[x y z] = size(im);
for ii = 2:nargin
    varargout{ii-1} = imresize(varargin{ii}, [x y]);
end




