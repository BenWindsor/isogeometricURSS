function dersv = perbspfunder( U, u, p, nders )
% Analogous to basisfunder.m but for periodic b splines and only up to
% degree 1
% INPUT
% U=knot vector
% u=eval points
% p=degree
% nders=number of derivatives

dersv = zeros(numel(u), nders+1, p+1);

% Iterate through eval points
for j=1:numel(u)
    % Find the non-zero periodic basis funs at u
    nonZeroBasis=numperbsp(U,u(j),p); 
    % Iterate through each basis function
    for k=1:(p+1)
        % Iterate through each derivative
        for l=1:(nders+1)
            if l==1
                dersv(j,l,k)=periodicBasisEval(U,u(j),nonZeroBasis(k),p);
            elseif l==2
                dersv(j,l,k)=periodicSingleBasisDerivEval(U,u(j),nonZeroBasis(k),p);
            else
                if p>1
                    dersv(j,l,k)=periodicSingleBasisSecondDerivEval(U, u(j), nonZeroBasis(k),p);
                elseif p==1
                    % If degree 1 just set second deriv to zero
                    dersv(j,l,k)=0;
                end
            end
        end
    end
end
        

end

