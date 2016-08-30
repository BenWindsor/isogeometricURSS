function val = zeroDegreeEval( U, u, i )
% Return the value of the spline N_{i,0}(u) on knot vector U

val=zeros(numel(u), 1);
for j=1:numel(u)
    if u(j)>=U(i) && u(j)<U(i+1)
        val(j)=1;
    else
        val(j)=0;
    end
end
end

