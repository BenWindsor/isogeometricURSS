function operator = localPeriodicOperator( U, u, p, i )

%One less element than knots?
elemNum=numel(U)-1; 

operator=zeros(p+1, p+1);

if k<=(p-1)
    operator(p+1, p+1)=1;
    for row=1:p
        for column=1:row
            operator(row,column)=alphaSubA(??);
        end
    end
elseif k>=p || k<=(elemNum-p+1)
    operator=diag(ones(p+1,1));
else
    operator(1, 1)=1;
    for row=2:(p+1)
        for column=row:p+1
            operator(row,column)=alphaSubA(??);
        end
    end
end


end

