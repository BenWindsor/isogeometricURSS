function perbspplot( pernurb, numPoints )
% plots a periodic basis spline
hold on;
func=@(x)(periodicSplineCurveEval(pernurb.knots,x,pernurb.order-1,pernurb.coefs));
points=func(linspace(pernurb.knots(1),pernurb.knots(end)-0.00001,numPoints));
scatter(points(1,:), points(2,:));

end

