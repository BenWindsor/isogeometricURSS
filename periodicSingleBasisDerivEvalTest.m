%Test periodic basis spline derivative
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
ctrl=[1 2 3 4 5];
crv=perbspmak(ctrl, U);
% 
% hold on;
% point=0.35;
% evalPoint=periodicSplineCurveEval(U, point, 2, ctrl);
% scatter(point, evalPoint);
% %plot the basis function in question
% perbspplot(crv,50);
% %plot the derivative tangent line to the point chosen above as u
% fplot(@(x)(periodicSplineCurveDerivEval(crv, point)*(x-point)+evalPoint), [0, 0.999]);
% 
% figure
% ctrl=[1 2 3 4 5; 0 1 -1 4 3];
% crv=perbspmak(ctrl, U);
% 
% hold on;
% perbspplot(crv, 50);
% 
% deriv=periodicSplineCurveDerivEval(crv,point);


% hold on;
% U=[0 0 0 0 0.2 0.4 0.6 0.8 1 1 1 1];
% elem=1;
% fplot(@(x)(periodicBasisEval(U, x, elem, 3)), [0 1]);
% point=0.19;
% eval=periodicBasisEval(U, point, elem, 3);
% deriv=periodicSingleBasisDerivEval(U, point, elem, 3);
% fplot(@(x)(deriv*(x-point)+eval), [0 1]);

hold on;
xHandle = @(x)(cos(2*pi*x));
yHandle = @(x)(sin(2*pi*x));
crv=periodicCurveInterpolate(15, 2, xHandle, yHandle);
perbspplot(crv, 100);
point = 0.0;
eval=perbspeval(crv, point);
tangent=@(x)(periodicSplineCurveDerivEval(crv, x)/norm(periodicSplineCurveDerivEval(crv, x)));
tangentEval=tangent(point);
quiver(eval(1), eval(2), tangentEval(1), tangentEval(2));
delta=0.001;
tangentDeriv=@(x)((tangent(x+delta)-tangent(x))/delta);
normal=@(x)(-tangentDeriv(x)/norm(tangentDeriv(x)));
normalEval=normal(point);
quiver(eval(1), eval(2), normalEval(1), normalEval(2));

