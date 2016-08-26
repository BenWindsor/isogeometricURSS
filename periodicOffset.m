function dist = periodicOffset( U, a, b )
% Returns the distance from a to b going periodically only to the right
% e.g if U=[0 1], a=0.9 and b=0.3 this should return 0.4
% INPUT
% U=knot vector
% a=first point
% b=second point

startPoint=U(1);
endPoint=U(end);

if a>endPoint || b>endPoint || a<startPoint || b<startPoint
    error('Point out of range');
end

if a<=b
    dist=b-a;
else
    dist=(endPoint-a)+(b-startPoint);
end


end

