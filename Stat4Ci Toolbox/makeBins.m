function bins = makeBins(data, nbCut)% Fr�d�ric Gosselin (frederic.gosselin@umontreal.ca), 09/03/2007temp = sort(data);cutoffs = zeros(1,nbCut+1);cutoffs(1) = temp(1)-1;for ii = 1:nbCut,	cutoffs(ii+1) = temp(ii*floor(length(temp)/nbCut));endtemp2 = data;for ii = 1:nbCut,	temp2(data>cutoffs(ii) & data<=cutoffs(ii+1)) = (ii-1) / (nbCut-1) - .5;endbins = temp2;