function varargout = overlay_cluster(varargin)

%overlay_cluster   Performs a cluster test and output a zscore map over an image. 
% 
%   out = overlay_cluster(im,zscores) Scales the zscore map to X% of
%   the total contrast of the image with a JET colormap. If zscores is a 3D
%   matrix of size x, y, z, the function tests z 2D matrices.  
% 
%   [out thresholded] =overlay_cluster(im,zscores,tC,k) This applies
%   the cluster test, threshold the zscore maps and scale the thresholded
%   zscores to X% of the image contrast with a JET colormap. If tC has more
%   than one element, it applies a two-tail cluster tests.
% 
% 
% 
% Nicolas Dupuis-Roy
% 2012-10-01
%
% See also overlay_pixel and find_sig_cluster

% if nargout>(nargin-2)
%     help overlay_cluster
%     error('wrong number of input/output arguments')
% end

% Assign input argument to a variable and put it in the right format
im = varargin{1}; im = stretch(double(im));
[x y z] = size(im);
if z==1
    im = repmat(im,[1 1 3]); %if not, make it 3D
end
zscores = varargin{2};

if nargin>2
    % Get input arguments
    tC = varargin{3}; 
    k = varargin{4}; 
    ndimtC = numel(tC);

    % Modify the k parameter to account for image size modifications
    [x2 y2 z2] = size(zscores);
    ndim = ndims(zscores);
    k_ = k^(1/ndim);
    k_new = (k_*x/x2) * (k_*y/y2);
    if ndim>2
        for ii = 3:ndims(zscores)
            k_new = k_new * k_;
        end
    end
    fprintf('Original k value: %d\n. Used k value: %d\n\n', k, k_new)
    zscores = sameSize(im, zscores); %Make it the same size as im (x- y- dimensions only)

    % Find the significant clusters
    sig_clust = find_sig_cluster(zscores,k_new,tC);
    for jj = 1:z2
        map{jj} = false(x,y,ndimtC);
        contourMap{jj} = false(x,y,ndimtC);
        c = 0;
        for ii = ndimtC:-1:1
            c = c + 1;
            map{jj}(:,:,c) = sig_clust(:,:,jj)==ii;
            contourMap{jj}(:,:,c) = edge(map{jj}(:,:,c),'canny')==0;
        end
    end
end
zscores = sameSize(im, zscores); %Make it the same size as im (x- y- dimensions only)

% Do the rest with the overlay function
%     - This function needs thresholded logical maps when nargout==2
if nargout==1
    varargout{1} = overlay(im,zscores);
    return;
elseif nargout>2
    varargout{1} = zeros(x,y,3,z2);
    varargout{2} = varargout{1};
    for jj = 1:z2
        [varargout{1}(:,:,:,jj) varargout{2}(:,:,:,jj)] = overlay(im,zscores(:,:,jj),map{jj},contourMap{jj});
%         [varargout{1}(:,:,:,jj) varargout{2}(:,:,:,jj)] = overlay(im,zscores(:,:,jj),map{jj});
    end
    if nargout==3
        varargout{3} = zeros(x,y,z2);
        for jj = 1:z2
            varargout{3}(:,:,jj) = map{jj};
        end
    end
end
