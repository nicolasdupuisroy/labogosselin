function EmailResults(Email, Title, Message, Files)
    setpref('Internet','E_mail', 'Viscog@UMontreal.CA');
    setpref('Internet','SMTP_Server', 'smtp.umontreal.ca');
    sendmail(Email, Title, Message, Files)
end
