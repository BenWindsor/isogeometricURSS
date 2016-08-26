function val = periodicBasisEval( U, u, elem, p )
% Returns the value of the i'th degree-p periodic basis spline on the knot vector U at
% point u
% INPUT:
% U = knot vector
% u = eval point
% elem=element which periodic spline eminates from
% p = degrree

% Which element u lies in
evalElem=findSpan(U,u)-p;
% Total elements
elemNum=numel(U)-2*p-1;

%WARNING: wont work for p not 2
if p==2
    % Case of no periodicity to account for
    if elem<=(elemNum-p)
        normalVals=localBasisSplineVectorEval(U, u, elem, p);
        operator=localPeriodicOperator(U, u, p, elem);
        newVals=operator*normalVals;
        val=newVals(3);
        % Case where 'end' functions overflow to 'beginning'
    else
        %The function overflowing into element 1
        if elem==(elemNum-1)
            if (elem+1)==evalElem
                normalVals=localBasisSplineVectorEval(U, u, evalElem, p);
                operator=localPeriodicOperator(U, u, p, evalElem);
                newVals=operator*normalVals;
                val=newVals(2);
            elseif evalElem==1
                normalVals=localBasisSplineVectorEval(U, u, evalElem, p);
                operator=localPeriodicOperator(U, u, p, evalElem);
                newVals=operator*normalVals;
                val=newVals(1);
            else
                fprintf('here\n');
                normalVals=localBasisSplineVectorEval(U, u, elem, p);
                operator=localPeriodicOperator(U, u, p, elem);
                newVals=operator*normalVals;
                val=newVals(3);
            end
            %The function overflowing into element 1 and 2
        elseif elem==(elemNum)
            if evalElem==1
                normalVals=localBasisSplineVectorEval(U, u, evalElem, p);
                operator=localPeriodicOperator(U, u, p, evalElem);
                newVals=operator*normalVals;
                val=newVals(2);
            elseif evalElem==2
                normalVals=localBasisSplineVectorEval(U, u, evalElem, p);
                operator=localPeriodicOperator(U, u, p, evalElem);
                newVals=operator*normalVals;
                val=newVals(1);
            else
                fprintf('here\n');
                normalVals=localBasisSplineVectorEval(U, u, elem, p);
                operator=localPeriodicOperator(U, u, p, elem);
                newVals=operator*normalVals;
                val=newVals(3);
            end
        end
    end
end

end
    




