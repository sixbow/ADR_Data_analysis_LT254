function [x_intersect,y_intersect] = findintersect_SdB(fh1,fh2,xguess)
% findintersect_SdB takes in 2 function handles of anonymous functions.
% you must provide a guess of x 
% Returns: x intersect near guess.
fming = @(x)fh1(x)-fh2(x);
x_intersect = fzero(fming,xguess);
y_intersect = fh1(x_intersect);
end

