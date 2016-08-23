function val = zeroDegreeEval( U, u, i )
% Return the value of the spline N_{i,0}(u) on knot vector U

if u>=U(i) && u<U(i+1)
    val=1;
else
    val=0;

end

