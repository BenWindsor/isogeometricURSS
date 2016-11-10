function val = periodicSplineCurveSecondDerivEval( perbsp, u)
% Returns the value of the second derivative of the spline curve at point u 
% INPUT:
% perbsp = periodic spline curve
% u = points to eval
% so each column is a coordinate point

U=perbsp.knots;
p=perbsp.order-1;
ctrl=perbsp.coefs;

% if p~=2
%     error('Can only handle degree 2');
% end

%If u is a cell convert to matrix
if iscell(u)
    u=cell2mat(u);
end

val=zeros(numel(ctrl(:,1)), numel(u));
elems=numel(U)-2*p-1;
for i=1:numel(u)
    for j=1:elems
        val(:,i) = val(:,i) + ctrl(:,j)*periodicSingleBasisSecondDerivEval(U, u(i), j, p);
    end
end

end

