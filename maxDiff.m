function diff = maxDiff( firstFunc, secondFunc, points )
% Compute the maximum difference between the first and second Function
% over the set of points given
% INPUT:
% firstFunc = function handle to first function
% secondFunc = function handle to second function
% points = points to test at

diffs=zeros(numel(points), 1);

for i=1:numel(points)
    diffs(i)=abs(firstFunc(points(i))-secondFunc(points(i)));
end

diff=max(diffs);

end

