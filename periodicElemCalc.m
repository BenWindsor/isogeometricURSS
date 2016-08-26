function newElem = periodicElemCalc( U, p, currentElem, additional )
% Determines which element you end up in if advancing +additional elements
% from currentElem, with periodic conditions in knot vector U
% INPUT
% U=knot vector
% p=degree
% currentElem=current element starting point
% additional=how many to jump forward

%Total element number
elemNum=numel(U)-2*p-1;

%IMPLEMENT: dealing with if additional is bigger than the number of
%elements e.g. adding 7 elements forward on a 5 elem space will currently
%still give overshoot

%If it stays in range
if (currentElem+additional)<=elemNum && (currentElem+additional)>=1
    newElem=currentElem+additional;
%If it goes off right hand side
elseif (currentElem+additional)>elemNum
    total=currentElem+additional;
    newElem=total-elemNum;
% If it goes off left hand side
else
    total=currentElem+additional;
    newElem=total+elemNum; %note total will be <=0 here
end

end

