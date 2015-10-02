function EmailError(email, err)
    setpref('Internet','E_mail', 'Masters@UMontreal.CA');
    setpref('Internet','SMTP_Server', 'smtp.umontreal.ca');
    [wtv, filename, ext] = fileparts(err.stack(1).file);
    sendmail(email, sprintf('ERROR in ''%s''', cat(2, filename, ext)), sprintf('%s\n\n%s', err.stack.file, err.getReport()), []);
end
