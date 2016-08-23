function pos = findSpan( U, u)
% Return index in knot vector that u falls into the interval of
% e.g. findSpan( [0 0.5 0.75 1], 0.53) = 2

for i=1:numel(U)
    if U(i)<=u && U(i+1)>u
        pos=i;
        break;
    end
end    


end

