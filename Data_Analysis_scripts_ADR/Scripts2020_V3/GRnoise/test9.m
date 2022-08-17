%Example input:
m = [1, 1; 2, 4; 3, 9];
filename = '' % do not write to file
t = matrix2latex(m, filename)
clipboard('copy',t) 