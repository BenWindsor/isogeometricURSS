function operator = localPeriodicOperatorDeg3(U, i)
% Returns the periodic operator on element i
% INPUT
% U=UNIFORM OPEN knot vector
% i=element index

p=3; %degree

%number accounts for repeated knots since is uniform open knot vector
elemNum=numel(U)-2*p-1; 

if i==1
    operator=[1/6 0 0 0; 2/3 2/3 0 0; 1/6 1/3 1 0; 0 0 0 1];
elseif i==2
    operator=[2/3 0 0 0; 1/3 1 0 0; 0 0 1 0; 0 0 0 1];
elseif i>=3 && i<=(elemNum-2)
    %Identity matrix
    operator=diag(ones(p+1,1));
elseif i==(elemNum-1)
    operator=[1 0 0 0; 0 1 0 0; 0 0 1 1/3; 0 0 0 2/3];
elseif i==elemNum
    operator=[1 0 0 0; 0 1 1/3 1/6; 0 0 2/3 2/3; 0 0 0 1/6];
end


end

