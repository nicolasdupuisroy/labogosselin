function varargout = ETA(varargin)

% ETA 
%   Start by initializing the function with: 
%         ETA('initialize', nbTrials)
%   Then, put the following line at the end of a loop: 
%         ETA(1) 
%   This prints the estimated time of arrival on the command window
%   based on a maximum sample of 30 iterations.  
% 
%   An alternative command line would be:
%         [year month day hour min sec] = ETA(0); 
%   This wont print the E.T.A. on the command window but will output the
%   same information.
% 
% Nicolas Dupuis-Roy
% 2012-10-01

persistent state estimated nbTrials

switch varargin{1}
    case 'initialize'
        nbTrials = varargin{2};
        state = 1;
        tic;
        estimated = [];
    otherwise
        print = varargin{1};
        state = state + 1;
        trial = state-1;
        t = toc;
        [year,month,day,hour,minute,second] = datevec(datenum(clock));
        
        % Estimate time of arrival (ETA)
        estimated = [estimated; t];
        if numel(estimated)>30
            estimated = estimated(end-29:end);
        end
        
        estimate = mean(estimated);
        sec_e = estimate*(nbTrials - trial);
        hour_e = floor(sec_e/3600);
        rem_sec = sec_e - hour_e*3600;
        min_e = floor(rem_sec/60);
        sec_e = round(rem_sec - min_e*60);

        sec_xtra = sec_e + second;
        min_e = min_e + floor(sec_xtra/60);
        sec_e = sec_xtra - 60*floor(sec_xtra/60);

        min_xtra = min_e + minute;
        hour_e = hour_e + floor(min_xtra/60);
        min_e = min_xtra - 60*floor(min_xtra/60);

        hour_xtra = hour_e + hour;
        day_e = day + floor(hour_xtra/24);
        hour_e = hour_xtra - 24*floor(hour_xtra/24);
        
        if print==1
            fprintf('ETA: %d-%d-%d %dh %dmin %dsec\n\n',year, month, day, hour_e, min_e, round(sec_e))
        end
        if nargout>0
            varargout{1} = year;
            varargout{2} = month;
            varargout{3} = day;
            varargout{4} = hour_e;
            varargout{5} = min_e;
            varargout{6} = round(sec_e);
        end
        tic
end

    
    
