function vals = localBasisSplineVectorEval( U, u, elem, p )
% Returns the vector N_elem(u) of basis functions with support on the knot
% span of that element, so elem=1 is the first non-empty knot section
% INPUT
% U=UNIFORM OPEN knot vector
% u=eval point
% elem=element of knot span
% p=degree

%NOTE: add error to handle if the point u isnt
%in the interval [U(i),U(i+1))

%IMPLEMENT: Do we really need elem provided? could use findSpan to get it
%instead??

vals=zeros(p+1, 1);
elemNum=numel(U)-2*p-1;
i=elem+p; %determine knot index

%Deal with left edge cases
if elem<=p
    for j=(i-p):i
        %Change index of vals to start from 1
        vals(j-(i-p)+1)=basisSplineEval(U,u,j,p);
    end
%Deal with right edge cases
elseif elem>=(elemNum-p)
    for j=(i-p):i
        %Change index of vals to start from 1
        vals(j-(i-p)+1)=basisSplineEval(U,u,j,p);
    end
%Normal middle case
else 
    start=i-p; 
    for j=1:numel(vals)
        vals(j)=basisSplineEval(U,u,start-1+j,p); %include the -1 since j starts at 1
    end
end

end

