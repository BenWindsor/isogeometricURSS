function dist = periodicOffset( U, a, b )
% Returns the distance from a to b going periodically only to the right
% e.g if U=[0 1], a=0.9 and b=0.3 this should return 0.4
% INPUT
% U=knot vector
% a=first point
% b=second point

startPoint=U(1);
endPoint=U(end);

if a>endPoint
    error('a too big');
end
if b>endPoint
    error('b too big');
end
if a<startPoint
    error('a too small');
end
if b<startPoint
    error('b too small');
end

if a<=b
    dist=b-a;
else
    dist=(endPoint-a)+(b-startPoint);
end


end

