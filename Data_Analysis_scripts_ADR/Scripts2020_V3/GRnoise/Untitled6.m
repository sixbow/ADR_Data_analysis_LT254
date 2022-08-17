LastName = {'Sanchez';'Johnson';'Li';'Diaz';'Brown'};             %                                        %
Weight = [176;163; 131; 133; 119];                                   %
T = table(Weight);                              %
T.Properties.RowNames = LastName;                                 %
table2latex(T);