function [r,p] = lsline3(varargin)

%lsline3 Plots the best-fitted line over the current figure. 
%   lsline2(x,y) Finds the best-fitted line over x and y and plots it over
%   the current figure. 
% 
%   lsline2(x,y,1) Finds the best-fitted line over x and y and plots it over
%   the current figure and prints the line equation plus RSQUARE on the
%   graph.
% 
% Nicolas Dupuis-Roy
% 2012-10-01
%
% Revision nov 2014 for MATLAB 2014a
% See also lsline


x = varargin{1};
y = varargin{2};

[r,p]=corrcoef(double(x(:)),double(y(:)));
stats = regstats(double(y(:)),double(x(:)));
xx = min(double(x(:))):.01:max(double(x(:)));
figure(gcf), hold on, plot(xx,xx*stats.beta(2) +  stats.beta(1),'--k')

if nargin==3
    Y = get(get(gcf,'children'), 'YLim');
    X = get(get(gcf,'children'), 'XLim');
    text(.7*range(X)+X(1),.8*range(Y)+Y(1),sprintf('y=%2.2fx+%2.2f\nR^2=%2.2f, p=%2.4f', stats.beta(2), stats.beta(1), r(1,2)^2, p(1,2)))
end
