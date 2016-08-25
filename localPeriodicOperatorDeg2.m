function operator = localPeriodicOperatorDeg2(U, i)
% Returns the periodic operator on element i
% INPUT
% U=UNIFORM OPEN knot vector
% i=element index e.g. working on [U(i), U(i+1)) knot span

%degree
p=2; 

%number accounts for repeated knots since is uniform open knot vector
elemNum=numel(U)-2*p-1; 

if i==1
    operator=[0.5 0 0; 0.5 1 0; 0 0 1];
elseif i>=2 && i<=(elemNum-1)
    %Identity matrix
    operator=diag(ones(3,1));
elseif i==elemNum
    operator=[1 0 0; 0 1 0.5; 0 0 0.5];
end


end
