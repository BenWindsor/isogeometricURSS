%Test derivatives 
%Plot a basis curve
U=[0,0,0, 0.2 0.4 0.6 0.8 1 1 1];
p=2;
hold on;
fplot(@(x)(basisSplineEval(U, x, 3, p)), [0, 0.9999]);

%Plot its derivative
fplot(@(x)(basisSplineDerivEval(U,x,3,p)), [0,0.9999]);

%Plot the tangent to that derivative i.e. second deriv at a point on
%derivative curve
point=0.3;
deriv=basisSplineSecondDerivEval(U, point, 3, p);
evalPoint=basisSplineDerivEval(U, point, 3, p);
fplot(@(x)(deriv*(x-point)+evalPoint), [0 0.9999]);

