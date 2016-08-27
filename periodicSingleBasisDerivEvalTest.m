%Test periodic basis spline derivative
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
ctrl=[1 2 3 4 5];
crv=perbspmak(ctrl, U);

hold on;
point=0.35;
evalPoint=periodicSplineCurveEval(U, point, 2, ctrl);
scatter(point, evalPoint);
%plot the basis function in question
perbspplot(crv,50);
%plot the derivative tangent line to the point chosen above as u
fplot(@(x)(periodicSplineCurveDerivEval(crv, point)*(x-point)+evalPoint), [0, 0.999]);

figure
ctrl=[1 2 3 4 5; 0 1 -1 4 3];
crv=perbspmak(ctrl, U);

hold on;
perbspplot(crv, 50);

deriv=periodicSplineCurveDerivEval(crv,point);