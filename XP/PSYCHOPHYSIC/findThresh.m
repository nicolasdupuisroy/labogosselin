function varargout = findThresh(varargin)

% 
% findThresh  This function finds the threshold of "input" at "crit"
%     according to a Weibull cumulative function. Note that output of the
%     Weibull function ranges from 0 to 1. Thus, if the input measurements
%     are not in this range, rescaling is needed (see third input argument).
%     When input is 2D, the second D is considered as the X values. 
% 
%     [thresh goodness model] = findThresh(input,crit) when
%     a third argument is inputed, the input is scaled according to the
%     min and max of value, and its returned thresholded value is rescaled
%     to the original domain.
% 
% 
%     [thresh goodness model] = findThresh(input,crit,[min max],flag) when
%     a fourth argument in inputed, a graphic is produced with fitting
%     data. Also, additionnal output arguments providing information on the
%     fitted model is given when requested.
% 
% See also FIT, FITTYPE
% 
% 
% Nicolas Dupuis-Roy, oct, 2012


% Prepare the input and assign variables.
input = varargin{1};
X = [1:length(input)]';
if size(input,1)==2;
    X = input(2,:);
    input = input(1,:);
    X = X(:);
    input = input(:);
elseif size(input2,2)==2;
    X = input(:,2);
    input = input(:,1);
end
input = input(:);
mnX = min(X);
mxX = max(X);
rX = range(X);
X = (X-mnX)/rX*4 + 1; %Rescaling of X so that it goes from 1 to 5

crit = varargin{2};
flag = false;
flag2 = false;
mnmx = [0 1];
if nargin>2
    mnmx = varargin{3};
    r = range(mnmx);
    input = (input-mnmx(1))/r; %First scale the input;
    flag2 = true;
    if nargin>3
        flag = true;
    end
end


% Fit the model
g = fittype('(1-exp(-(x/h).^b))');

[model goodness] = fit(X,input,g,'Startpoint', [7 3]);

% The model
x = linspace(0,numel(input)+10,10000);
y = model(x)*r+mnmx(1);


% Prepare and assign the output
tmp = x(find(y>=crit, 1,'first'));
varargout{1} = ((tmp - 1)/4)*rX + mnX; %Rescaling of X
if isempty(varargout{1})
    disp('Threshold out of range!')
    return
end
if flag2
    input = input*r+mnmx(1);
end
if nargout>1
    varargout{2} = goodness;
    if nargout>2
        varargout{3} = model;
    end
end


% Produce a graphic if requested
X = (X-1)/4*rX + mnX;
if flag
    % Prepare lines that points toward the threshold
    xxx1 = ones(1,100)*varargout{1};
    yyy1 = linspace(min(input),crit,100);
    xxx2 = linspace(min(X),varargout{1},100);
    yyy2 = ones(1,100)*crit;
    figure, plot(xxx1,yyy1,'k--'), hold on, 
    plot(xxx2,yyy2,'k--'), hold on, 
    
    x = linspace(1,5,1000);
    y = model(x)*r+mnmx(1);
    x = linspace(min(X),max(X),1000);

    % The actual data and fitted data
    plot(X,input,'.'), hold on, 
    plot(x,y,'-k')
    
    % Model's fit info
    txt1 = sprintf('thresh(%2.2f) = %2.4f',crit,varargout{1});
    text(x(1)+range(x)*.2, max(y)-.1*range(y), txt1);
    txt2 = sprintf('R^2 = %2.2f',goodness.rsquare);    
    text(x(end)-range(x)*.2, min(y)+.1*range(y), txt2);
    axis tight
end
