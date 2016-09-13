function val = periodicBasisEval( U, u, elem, p )
% Returns the value of the i'th degree-p periodic basis spline on the knot vector U at
% point u
% INPUT:
% U = knot vector
% u = eval point/s
% elem=element which periodic spline eminates from
% p = degrree

%IMPLEMENT: attempting to evaluate at last point on knot will return the
%value at first point of knot as it is periodic, allowing evaluation at
%e.g. 1 on a knot span [0 1] 

%convert cell of points to normal array if needed
if iscell(u)
    u=cell2mat(u);
end

% Total elements
elemNum=numel(U)-2*p-1;

val=zeros(numel(u),1);
%WARNING: wont work for p not 2
for j=1:numel(u)
    % Which element u lies in, evaluating at 0 if it is on the end knot
    if u(j)<U(end)
        evalElem=findSpan(U,u(j))-p;
    elseif u(j)==U(end)
        val(j)=periodicBasisEval(U,0,elem,p);
        break;
    end
    
    if p==1
        % Case of no periodicity to account for
        if elem<=(elemNum-p)
            normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
            operator=localPeriodicOperator(U, u(j), p, elem);
            newVals=operator*normalVals;
            val(j)=newVals(2);
            % Case where 'end' functions overflow to 'beginning'
        else
            if elem==elemNum
                if evalElem==1
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(1);
                else
                    normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
                    operator=localPeriodicOperator(U, u(j), p, elem);
                    newVals=operator*normalVals;
                    val(j)=newVals(2);
                    
                end
            end
        end
        
    elseif p==2
        % Case of no periodicity to account for
        if elem<=(elemNum-p)
            normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
            operator=localPeriodicOperator(U, u(j), p, elem);
            newVals=operator*normalVals;
            val(j)=newVals(3);
        % Case where 'end' functions overflow to 'beginning'
        else
            %The function overflowing into element 1
            if elem==(elemNum-1)
                if (elem+1)==evalElem
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(2);
                elseif evalElem==1
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(1);
                else
                    normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
                    operator=localPeriodicOperator(U, u(j), p, elem);
                    newVals=operator*normalVals;
                    val(j)=newVals(3);
                end
                %The function overflowing into element 1 and 2
            elseif elem==(elemNum)
                if evalElem==1
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(2);
                elseif evalElem==2
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(1);
                else
                    normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
                    operator=localPeriodicOperator(U, u(j), p, elem);
                    newVals=operator*normalVals;
                    val(j)=newVals(3);
                end
            end
        end
        
    elseif p==3
           % Case of no periodicity to account for
        if elem<=(elemNum-p)
            normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
            operator=localPeriodicOperator(U, u(j), p, elem);
            newVals=operator*normalVals;
            val(j)=newVals(4);
            % Case where 'end' functions overflow to 'beginning'
        else
            %The function overflowing into element 1
            if elem==(elemNum-2)
                if (elem+1)==evalElem
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(3);
                elseif (elem+2)==evalElem
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(2);
                elseif evalElem==1
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(1);
                else
                    normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
                    operator=localPeriodicOperator(U, u(j), p, elem);
                    newVals=operator*normalVals;
                    val(j)=newVals(4);
                end
            %The function overflowing into element 1 and 2
            elseif elem==(elemNum-1)
                if (elem+1)==evalElem
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(3);
                elseif evalElem==1
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(2);
                elseif evalElem==2
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(1);
                else
                    normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
                    operator=localPeriodicOperator(U, u(j), p, elem);
                    newVals=operator*normalVals;
                    val(j)=newVals(4);
                end
            %The function overflowing into element 1, 2 and 3
            elseif elem==(elemNum)
                if evalElem==1
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(3);
                elseif evalElem==2
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(2);
                elseif evalElem==3
                    normalVals=localBasisSplineVectorEval(U, u(j), evalElem, p);
                    operator=localPeriodicOperator(U, u(j), p, evalElem);
                    newVals=operator*normalVals;
                    val(j)=newVals(1);
                else
                    normalVals=localBasisSplineVectorEval(U, u(j), elem, p);
                    operator=localPeriodicOperator(U, u(j), p, elem);
                    newVals=operator*normalVals;
                    val(j)=newVals(4);
                end
            end
        end
    else
        error('only available for degree p=2 or p=3');
    end
end
end
    




