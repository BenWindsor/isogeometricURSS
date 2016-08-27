function pos = findSpan( U, u)
% Return index in knot vector that u falls into the interval of
% e.g. findSpan( [0 0.5 0.75 1], 0.53) = 2
%INPUT
% U=knot vector
% u=value or values to find knots for 

% NOTE: this is the same as the NURB Toolbox findspan.m but just starting
% counting from 1 not 0.

pos=zeros(size(u));
for j=1:numel(u)
    for i=1:numel(U)
        if U(i)<=u(j) && U(i+1)>u(j)
            pos(j)=i;
            break;
        end
    end
end

end

