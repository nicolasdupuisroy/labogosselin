function [Ci, nbIncluded] = buildCi2(cid, predVar, incCriteria)% BUILDCI Performs least-square multiple linear regression on the sampling noise % 	(explanatory variable) rebuilt using the function contained in the noise input % 	string and saved in a temporary M-file named MAKETHENOISE; and the elements of% 	the data input matrix (i.e., trials) satisfying two conditions (predictive % 	variable): (1) belonging to a column identified in the predVar vector (e.g., % 	predVar = [8, 10]) and (2) verifying the Boolean inclusion criterias specified % 	in the cell of strings incCriteria. In incCriteria, the expression 'data(n,:)' % 	refers to the nth data column.% 		E.g., [data, origins, info, noise, dataLabels] = readCid('AGBTSD04.cid');% 			  incCriteria{1} = '(data(2,:)==1)'; %2nd data columns = 1% 			  incCriteria{2} = '(data(2,:)==1)&(data(3,:)==1)'; %2nd and 3rd data columns = 1%	 		  predVar = [8, 10];% 			  [Ci, nbIncluded] = buildCi(data, noise, predVar, incCriteria);% %	BUILDCI goes through all of the data input matrix. Therefore, you should only pass to % 	this function the portion of the data input matrix that satisfies to the union of all% 	inclusion criteria. Both NOISE and data are outputs of the READCID function. If predVar % 	has m elements and incCriteria contains n elements, m * n classification images are computed. % 	The Ci output matrix contains them all; it has a size of d_1 * ... * d_z * m * n, where % 	d_x give the number of elements in dimension x, z is the number of dimensions of the % 	sampling noise, m is the number of elements in the predVar vector and n the number of % 	strings in the incCriteria input cell. The number of trials used to build Ci is outputed % 	in the scalar nbIncluded. The elements of the Ci are Z-scores.% %	Importantly, it is assumed that the sampling noise is uncorrelated between trials. %	This is true if the sampling noise is random. Least-square multiple regression then %	becomes Ci = ky'X, where Ci is a vectorized classification image, k is a scalar, y %	is a predictive variable vector, and X is a explanatory variable matrix. Up to the % 	k factor, this amounts to summing all sampling noise matrices multiplied by their % 	predictive variable values.% % See also READCID, DEMO4CI% % The Stat4Ci toolbox is free (http://mapageweb.umontreal.ca/gosselif/stat4ci.html); if you use % it in your research, please, cite us:%	Chauvin, A., Worsley, K. J., Schyns, P. G., Arguin, M. & Gosselin, F. (2004).  A sensitive %	statistical test for smooth classification images.% % Alan Chauvin & Fr�d�ric Gosselin (frederic.gosselin@umontreal.ca), 27/08/2004% Modified by Fr�d�ric Gosselin, 16/02/2006% to do: Test the inclusion algorithm% to do: Allow the use of the column names instead of 'C1' or 'data(1,:)'% for the incCriteria + regressor (transform into a cell?)% to do: test of z-scoring with noise + something already analyzed% to do: rewrite the help section (requires bwlabel from the image processing toolbox)noise = cid.noise;constants = cid.constants;eval(constants);data = cid.data;dataLabels = cid.dataLabels;[nbColumns, nbTrials] = size(data);nbPredVar = length(predVar);[nothing, nbIncCriteria] = size(incCriteria);id = fopen('makeTheNoise.m', 'w');fprintf(id, '%s', noise);fclose(id);criteriaMet = zeros(nbIncCriteria, nbTrials);for whichIncCriteria = 1:nbIncCriteria,	outputStr = incCriteria{whichIncCriteria};	for ii = nbColumns:-1:1,		outputStr = replaceStr(outputStr, sprintf('C%d',ii), sprintf('data(%d,:)',ii));		outputStr = replaceStr(outputStr, dataLabels{ii}, sprintf('data(%d,:)',ii));	end	criteriaMet(whichIncCriteria, :) = eval(outputStr);end%%%%%%%%%%%%%%%%%%%%%%%collapsed = ge(sum(eq(data, -99), 1), 1);index = ge(sum(criteriaMet, 1), 1);newIndex = criteriaMet;lIndex = bwlabel(index);nbSeg = max(lIndex);for ii = 1:nbSeg,    temp = eq(lIndex, ii);    toto = find(temp);	ind = toto(1);	    while collapsed(ind) == 1        if ind == 1            error('No seed provided.')        end        ind = ind - 1;        newIndex(ind) = 1;        endend% index = newIndex;index = ge(sum(newIndex, 1), 1);data = data(:, eq(index, 1));nbTrials = size(data, 2);criteriaMet = criteriaMet(:, eq(index, 1));% % % % % % % % % % % % % % % % % % % % % % % % % % %% fucks when incCriteria does not include everything% % % % % % % % % % % % % % % % % % % % % % % % % % %nbIncluded = zeros(1, nbIncCriteria);whichTrial = 1;oneTrial = data(:,whichTrial);theNoise = makeTheNoise(oneTrial);nbDim = length(size(theNoise));Ci = zeros([size(theNoise), nbPredVar, nbIncCriteria]);% should use a cell to identify automatically all the dimensions to include in one, independent Ciglobal indCiif (isempty(indCi) | indCi == 0)	nbCi = 1;	indCi = 0;else	nbCi = size(theNoise, indCi);end% should use a cell to identify automatically all the dimensions to include in one, independent CigrandMeans = zeros(nbCi, nbPredVar, nbIncCriteria);grandVars = zeros(nbCi, nbPredVar, nbIncCriteria);temp = 'Ci(';for ii = 1:nbDim, temp = [temp,':,'];, endtemp = [temp,'whichPredVar,whichIncCriteria)'];const = sprintf('%s = %s + oneTrial(predVar(whichPredVar)) * theNoise;', temp, temp);fprintf('\n%d%%', round(100*whichTrial/nbTrials));for whichIncCriteria = 1:nbIncCriteria,	if(criteriaMet(whichIncCriteria, whichTrial))		nbIncluded(whichIncCriteria) = nbIncluded(whichIncCriteria) + 1;		for whichPredVar = 1:nbPredVar,			eval(const)			for whichCi = 1:nbCi,				temp = 'theNoise(';				for ii = 1:nbDim,					if ii == indCi,						temp = [temp,sprintf('%d,', whichCi)];						else						temp = [temp,':,'];					end				end				temp = [temp(1:end-1),')'];				tempNoise = eval(temp);				grandMeans(whichCi, whichPredVar, whichIncCriteria) = grandMeans(whichCi, whichPredVar, whichIncCriteria) + oneTrial(predVar(whichPredVar)) * mean(tempNoise(:));				grandVars(whichCi, whichPredVar, whichIncCriteria) = grandVars(whichCi, whichPredVar, whichIncCriteria) + (oneTrial(predVar(whichPredVar)) * std(tempNoise(:)))^2;			end		end	endendfor whichTrial = 2:nbTrials,	back = '\b';	nbBack = length(num2str(round(100*(whichTrial-1)/nbTrials)));	for ii = 1:nbBack, back = [back, '\b']; end	eval(sprintf('fprintf(''%s%d%s'');', back, round(100*whichTrial/nbTrials),'%%'));	oneTrial = data(:,whichTrial);	theNoise = makeTheNoise(oneTrial);	for whichIncCriteria = 1:nbIncCriteria,		if(criteriaMet(whichIncCriteria, whichTrial))			nbIncluded(whichIncCriteria) = nbIncluded(whichIncCriteria) + 1;			for whichPredVar = 1:nbPredVar,				eval(const)						for whichCi = 1:nbCi,					temp = 'theNoise(';					for ii = 1:nbDim,						if ii == indCi,							temp = [temp,sprintf('%d,', whichCi)];							else							temp = [temp,':,'];						end					end					temp = [temp(1:end-1),')'];					tempNoise = eval(temp);					grandMeans(whichCi, whichPredVar, whichIncCriteria) = grandMeans(whichCi, whichPredVar, whichIncCriteria) + oneTrial(predVar(whichPredVar)) * mean(tempNoise(:));					grandVars(whichCi, whichPredVar, whichIncCriteria) = grandVars(whichCi, whichPredVar, whichIncCriteria) + (oneTrial(predVar(whichPredVar)) * std(tempNoise(:)))^2;				end			end		end	endendfor whichIncCriteria = 1:nbIncCriteria,	if(criteriaMet(whichIncCriteria, whichTrial))		for whichPredVar = 1:nbPredVar,						for whichCi = 1:nbCi,				temp = 'Ci(';				for ii = 1:nbDim,					if ii == indCi,						temp = [temp,sprintf('%d,', whichCi)];						else						temp = [temp,':,'];					end				end				temp = [temp,'whichPredVar,whichIncCriteria)'];				temp = [temp, '=(', temp, '- grandMeans(whichCi,whichPredVar,whichIncCriteria)) / sqrt(grandVars(whichCi,whichPredVar,whichIncCriteria));'];				eval(temp)			end		end	endendCi = squeeze(Ci);%---------------------------------------------------------------%---------------------------------------------------------------function outputStr = replaceStr(inputStr, findThis, changeFor)index = findstr(findThis, inputStr);NN = length(index);m1 = length(findThis);m2 = length(changeFor);inSize = length(inputStr);index(NN+1) = inSize+1;if(index(1)~=1)	outputStr(1:index(1)-1) = inputStr(1:index(1)-1);endfor ii = 1:NN,	jj1 = index(ii);	jj2 = index(ii+1);	temp = inputStr(jj1+m1:jj2-1);	pp = length(temp);	kk = (ii-1) * m2 + (jj1 - (ii - 1) * m1);	outputStr(kk:kk+m2-1) = changeFor;	outputStr(kk+m2:kk+m2+pp-1) = temp;end