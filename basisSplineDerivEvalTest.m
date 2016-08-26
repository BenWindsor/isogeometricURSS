%Test derivatives 
U=[0,0,0, 0.2 0.4 0.6 0.8 1 1 1];
p=2;
hold on;
for i=1:7
    fplot(@(x)(basisSplineEval(U,x,i,p)), [0,0.9999]);
end

point=0.64;
deriv=basisSplineDerivEval(U, point, 5, p);
evalPoint=basisSplineEval(U, point, 5, p);

fplot(@(x)(deriv*(x-point)+evalPoint), [0 0.9999]);

