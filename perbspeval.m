function val = perbspeval( perbsp, u )
% evaluate the perbsp structure at point u
func=@(x)(periodicSplineCurveEval(perbsp.knots,x,perbsp.order-1,perbsp.coefs));


val=func(u);

end

