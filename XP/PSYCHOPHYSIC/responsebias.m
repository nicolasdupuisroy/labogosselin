function c = responsebias(hit_rate,false_alarm)
% function c = responsebias(hit_rate, false_alarm)
%
% Compute c given hit_rate and false_alarm
% 2012-10-01 Nicolas Dupuis-Roy

if (hit_rate == 1)
	hit_rate = .999;
elseif (hit_rate == 0)
	hit_rate = .001;
end

if (false_alarm == 1)
	false_alarm = .999;
elseif (false_alarm == 0)
	false_alarm = .001;
end


c = -(invNormCum(hit_rate,0,1) + invNormCum(false_alarm,0,1))/2;

