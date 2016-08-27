function B = numperbsp( U, u, p )
% Analogue of numbasisfun.m from NURBS toolbox, it tells you the numbers of
% the starting elements for non-zero basis functions at point u in knot
% vector U of degree p 

B=zeros(numel(u), p+1);
numElem=numel(U)-2*p-1; %Number of elements in span

for i=1:numel(u)
    
    %Find what element the point u is in
    span=findSpan(U,u(i));
    if span<=p
        elem=1;
    elseif span>=(numel(U)-p)
        elem=numElem;
    else 
        elem=span-p;
    end
    
    
    for j=1:(p+1)
        B(i,j)=periodicElemCalc(U,p,elem,-(j-1));
    end

end

