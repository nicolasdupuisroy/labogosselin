function out = find_sig_cluster(zscores,k, tC)

%find_sig_cluster    Find significant clusters in a 2D or 3D zscore map.
%   out = find_sig_cluster(zscores,k,tC) find clusters within the thresholded 
%   zscores map that are larger or equal to k. If tC has two values, a two tail
%   cluster test is applied. In this case, the significant clusters are
%   labeled 1 and 2 refering to negative and positive tC. If tC
%   has more than one element, it does a two-tail pixel test. 
% 
% Nicolas Dupuis-Roy
% 2012-10-01
%
%   See also bub_overlay_pixel and bub_overlay_cluster

ndim = numel(size(zscores)>1);
tC = sort(tC);
ndimTc = numel(tC);
thresh = cell(1,ndimTc);
label_tag = 1:ndimTc;

out = zeros(size(zscores), 'uint8');
thresh{1,end} = zscores>=tC(end);
if ndimTc>1
    thresh{1,1} = zscores<=tC(1);
end

if ndim==2
    conn = 4;
elseif ndim==3
    conn = 6;
end

for ii = 1:ndimTc
    [labels, nlabels] = bwlabeln(thresh{1,ii},conn);

    for jj = 1:nlabels
        S = sum(labels(:)==jj);
        if S>=k
            fprintf('Significant cluster of size %d found (k=%d)\n', S, k)
            out = uint8(labels==jj)*label_tag(ii) + out;
        else
            fprintf('Cluster of size %d found (k=%d)\n', S, k)
        end
    end
end
