function operator = localPeriodicOperator( U, u, p, i )
% Return local operator that transforms a b-spline into a periodic spline
% INPUT
% U=UNIFORM OPEN knot vector
% u=eval point or knot?? not sure yet
% p=degree
% i=element index e.g. [U(i), U(i+1)) interval


%Accounts for repeated knots since open vector
elemNum=numel(U)-2*p-1; 



%TEMPORARY IMPLEMENTATION
if p==2
    operator=localPeriodicOperatorDeg2(U, i);
elseif p==3
    operator=localPeriodicOperatorDeg3(U, i);
end

%FULL IMPLEMENTATION
%operator=zeros(p+1, p+1);
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

