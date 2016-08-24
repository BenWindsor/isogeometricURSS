function alpha = alphaSubA(U, u, p, A)
% Return the value of \alpha_A from the bezier extraction Hughes paper pg4
% equation(10)
% INPUT
% U=knot vector (size n+p+1 in paper notation)
% u=new knot falling in [U(k), U(K+1)) interval to be found by findSpan
% p=degree
% A=subscript of alpha

k=findSpan(U, u);

if 1<=A && A<= (k-p)
    alpha=1;
elseif (k-p-1)<=A && A<=k
    alpha=(u-U(A))/(U(A+p)-U(A));
elseif A>=(k+1)
    alpha=0;
end

end

