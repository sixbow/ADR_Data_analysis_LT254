% This script can be run after S21analysis_GRV5.m 

% Saves the resonator Qi, Ql, Qc, Fres, P_int, P_read,
% KIDnumber, and the errrors of  Qi, Qc, Ql and Fres in a format that can
% easily be imported using Python (numpy)

output.KIDnumber = [KID.KIDnumber];
output.ReadPower = [KID.ReadPower];
output.InternalPower = [KID.InternalPower];
output.Fres = [KID.Fres];
output.Ql = [KID.Ql];
output.Qi = [KID.Qi];
output.Qc = [KID.Qc];

j = 1;
for i = 1:length([KID.KIDnumber])
    errors = [KID.Paramerrors];
    output.Fres_error(i) = errors(j);
    output.Ql_errror(i) = errors(j+2);
    output.Qi_errror(i) = errors(j+3);
    output.Qc_errror(i) = errors(j+4);
    
    j = j + 5;
end

save('aSi_250C_lowP','output')