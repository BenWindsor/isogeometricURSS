function vals = localNurbVectorEval( U, u, i, p, W )
% Returns the vector N_i(u) of basis functions with support on the knot
% span [U(i), U(i+1)). This will be a vector of length P+1.
% INPUT
% U=knot vector
% u=eval point
% i=knot index for start of span
% p=degree
% W=knot vector

%NOTE: add error to handle if the point u isnt
%in the interval [U(i),U(i+1))

%NOTE: currently i is not used, should correct

vals=zeros(p+1, 1);
start=i-p;

for j=1:numel(vals)
    vals(j)=nurbBasisEval(U,u,start-1+j,p,W); %include the -1 since j starts at 1
end

end

