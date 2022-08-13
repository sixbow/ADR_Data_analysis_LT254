function gensec(tabs,title)
%gensec makes the comments lines for seperation section in your code.
tabstr= [''];
for i=1:tabs
    tabstr = [tabstr '\t'];
end
dashstr= [''];
dash = 40+22;
for i=1:(dash-4*tabs-length(char(title)))
    dashstr = [dashstr '-'];
end
formatstr_plus = [ append(tabstr) '%s+++++|Begin:%s' append(dashstr) '\n'];
begin_str = sprintf(formatstr_plus,'%',title);
formatstr_plus = [ append(tabstr) '%s-------/End:%s' append(dashstr) '\n'];
end_str = sprintf(formatstr_plus,'%',title);
tot_str = append(begin_str,end_str);
clipboard('copy',tot_str)  
end

