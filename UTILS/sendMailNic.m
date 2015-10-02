function sendMailNic(varargin)

if nargin==1
    error('this function needs at least 2 arguments')
end
message = [];
if nargin>1
    email = varargin{1};
    goodOrBad = varargin{2};
end
if nargin==3
    message = varargin{3};
end
    
% Prepare    
mail = 'matlabdebug@gmail.com';
password = 'jouissance';
setpref('Internet','SMTP_Server','smtp.gmail.com')
setpref('Internet','SMTP_Username',mail)
setpref('Internet','SMTP_Password',password)

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Send message to varargin{1}
if goodOrBad==1
    title = 'Eh oh';
    if isempty(message)
        message = 'Ton analyse est terminé mon coco!';
    end
else
    title = 'Bad news';
    if isempty(message)
        message = 'Something bad happened';
    end
end
sendmail(email,title,message)