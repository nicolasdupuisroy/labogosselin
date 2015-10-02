function visualAngle

% visualAngle     This function implement the formula to compute the visual
% angle (deg). Enter two values, leave one blank and the answer will be
% inputed in the dialog. If all three values are inputed, than the value
% for the visual angle could be validated.

name='Visual Angle Formula (Leave blank if unknown)';
prompt={'Visual angle (in degrees):','Size of object (in cm):', 'Distance to object (in cm):'};
numlines=1;
defaultanswer={'1.03','2.3',''};
answer=inputdlg(prompt,name,numlines,defaultanswer);

empty = false(1,3);
for ii = 1:3
    answer{ii} = str2num(answer{ii});
    empty(ii) = isempty(answer{ii});
end

if sum(empty)>1
    error('You have more than one unknown parameters (or you entered a string)!')
end
if sum(empty)==0
    ButtonName = questdlg('Seems that you know all the answer. Would you like me to check if your computation is correct?', ...
                         'Validation', ...
                         'Yes', 'No', 'Yes');
   switch ButtonName,
     case 'No',
      return;
     case 'Yes',
        ButtonName2 = questdlg('Which value would you like to validate?', ...
                         'Validation', ...
                         'Angle', 'Size', 'Distance', 'Angle');
        switch ButtonName2
            case 'Angle'
                err = answer{1} - 2*atan(answer{2}/2/answer{3})*180/pi;
                err = abs(err/answer{1}*100);
                answer{1} = 2*atan(answer{2}/2/answer{3})*180/pi;
            case 'Size'
                err = answer{2} - 2*tan(answer{1}/180*pi/2)*answer{3};
                err = abs(err/answer{2}*100);
                answer{2} = 2*tan(answer{1}/180*pi/2)*answer{3};
            case 'Distance'
                err = answer{3} - answer{2}/2/tan(answer{1}/180*pi/2);
                err = abs(err/answer{3}*100);
                answer{3} = answer{2}/2/tan(answer{1}/180*pi/2);
        end
        questdlg(sprintf('The error on the %s is of %3.2f%%. The %s will be corrected on the next input dialog!',...
            ButtonName2,err,ButtonName2),...
             'Validataion Output',...
             'Ok', 'Ok');
   end
end
if find(empty)==1
    answer{1} = 2*atan(answer{2}/2/answer{3})*180/pi;
elseif find(empty)==2
    answer{2} = 2*tan(answer{1}/180*pi/2)*answer{3};
elseif find(empty)==3
    answer{3} = answer{2}/2/tan(answer{1}/180*pi/2);
end

for ii = 1:3
    defaultanswer{ii}=num2str(answer{ii});
end
inputdlg(prompt,name,numlines,defaultanswer);
