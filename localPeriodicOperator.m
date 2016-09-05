function operator = localPeriodicOperator( U, u, p, elem )
% Return local operator that transforms a b-spline into a periodic spline
% INPUT
% U=UNIFORM OPEN knot vector
% u=eval point or knot?? not sure yet
% p=degree
% elem=element index 





% TEMPORARY IMPLEMENTATION
if p==1
    operator=[1 0; 0 1];
elseif p==2
    operator=localPeriodicOperatorDeg2(U, elem);
elseif p==3
    operator=localPeriodicOperatorDeg3(U, elem);
end

% FULL IMPLEMENTATION
% Accounts for repeated knots since open vector
% elemNum=numel(U)-2*p-1; 
% operator=zeros(p+1, p+1);
% if i<=(p-1)
%     operator(p+1, p+1)=1;
%     for row=1:p
%         for column=1:row
%             operator(row,column)=alphaSubA(??);
%         end
%     end
% elseif i>=p && i<=(elemNum-p+1)
%     operator=diag(ones(p+1,1));
% else
%     operator(1, 1)=1;
%     for row=2:(p+1)
%         for column=row:p+1
%             operator(row,column)=alphaSubA(??);
%         end
%     end
% end


end

