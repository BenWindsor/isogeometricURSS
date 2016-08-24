function invWMat = invWeightMatrix( W )
% Takes a set of weights and returns diag(W)^{-1}

invW=zeros(numel(W),1);

for i=1:numel(W)
    invW(i)=1/(W(i));
end

invWMat=diag(invW);


end

