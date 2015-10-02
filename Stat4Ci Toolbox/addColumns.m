function cid = addColumns(cid, cols, names)% ADDCOLUMNS Adds cols, a M x N matrix, M being the number of columns added and N the number of % lines in cid.data, at the end of the cid.data matrix and adds names, a M element cell, at the % end of the cell cid.dataLabels.  Both modifications are required to use the makeCid function.% % See also READCID, DEMO4CI% % The Stat4Ci toolbox is free (http://mapageweb.umontreal.ca/gosselif/stat4ci.html); if you use % it in your research, please, cite us:%	Chauvin, A., Worsley, K. J., Schyns, P. G., Arguin, M. & Gosselin, F. (2004).  A sensitive %	statistical test for smooth classification images.% % Fr�d�ric Gosselin (frederic.gosselin@umontreal.ca), 09/03/2007cid.data(end+1:end+size(cols, 1), :) = cols;cid.dataLabels = [cid.dataLabels, names];