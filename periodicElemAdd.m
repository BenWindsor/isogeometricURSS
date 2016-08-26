function newElem = periodicElemAdd( U, p, currentElem, additional )
% Determines which element you end up in if advancing +additional elements
% from currentElem, with periodic conditions in knot vector U
% INPUT
% U=knot vector
% p=degree
% currentElem=current element starting point
% additional=how many to jump forward

%Total element number
elemNum=numel(U)-2*p-1;
if (currentElem+additional)<=elemNum
    newElem=currentElem+additional;
else
    total=currentElem+additional;
    newElem=total-elemNum
end

end

