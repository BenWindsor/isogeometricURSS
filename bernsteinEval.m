function val = bernsteinEval( a, p, u )
% Returns the value of the Bernstein polynomial B_{a,p}(u) as detailed
% in the isogeometric finite element bezier extraction paper 

% NOTE: the parameter p seems to vary in literature, here we have that p=1
% is the first bernstein polynomial, constant y=1. Sometimes this is
% denoted with p=0. The wikipedia uses p=0 but the Hughes paper named above
% uses p=1 so I am using that convention for now as that is what they use
% to implement their bezier extraction operator.
if a==1 && p==1
    val = 1;
elseif a<1 | a>p+1
    val = 0;
else
    val = (1-u)*bernsteinEval(a, p-1, u) + u*bernsteinEval(a-1, p-1, u);
end

