function S21 = totalS21(Filter)
% filter is a n array of cells with rw, ncols with S21 data of each step.
% FUnction calcualtes the total S21
S21 = Filter{1}(1,:);
for n = 2 : length(Filter)
    S21 = S21 .* Filter{n}(1,:);
end

end