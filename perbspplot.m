function perbspplot( pernurb, numPoints )
% plots a periodic basis spline
hold on;
func=@(x)(periodicSplineCurveEval(pernurb.knots,x,pernurb.order-1,pernurb.coefs));



if pernurb.dim==1
    fplot(func, [pernurb.knots(1) pernurb.knots(end)-0.0001]);
elseif pernurb.dim==2
    points=func(linspace(pernurb.knots(1),pernurb.knots(end)-0.00001,numPoints));
    scatter(points(1,:), points(2,:));
else
    error('Can only handle dim 1 or 2');
end

end

