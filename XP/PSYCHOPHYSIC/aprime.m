function ap = aprime(hit, false_alarm)
% function ap = aprime(hit, false_alarm)
%
% Computes A' given hit and false_alarm
% Formulas from Snoodgrass and Corwin (1988). J. Exp. Psychol., 117, p.38
% 26/04/2009 Frederic Gosselin

if hit>=false_alarm,
   ap = .5+((hit-false_alarm)*(1+hit-false_alarm))/(4*hit*(1-false_alarm));
else
   ap = .5-((false_alarm-hit)*(1+false_alarm-hit))/(4*false_alarm*(1-hit));
end
