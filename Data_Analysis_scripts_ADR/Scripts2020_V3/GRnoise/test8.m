LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};             %
A = [1 1;1 1;1 1;1 1;1 1]
T = table(A);                              %
A.Properties.RowNames = LastName;                                 %
table2latex(A);