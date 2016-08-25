function vals = localBasisSplineVectorEval( U, u, i, p )
% Returns the vector N_i(u) of basis functions with support on the knot
% span [U(i), U(i+1)). This will be a vector of length P+1.
% INPUT
% U=UNIFORM OPEN knot vector
% u=eval point
% i=knot index for start of span 
% p=degree

%NOTE: add error to handle if the point u isnt
%in the interval [U(i),U(i+1))

vals=zeros(p+1, 1);

% %Give correct start, taking into account uniform knot vectors
% if i>=1 && i<=p
%     %So the 'first' knot to eval is the first non-empty interval
%     start=1;
% else
%     start=i-p;
% end

start=i-p;

for j=1:numel(vals)
    vals(j)=basisSplineEval(U,u,start-1+j,p); %include the -1 since j starts at 1
end

end

