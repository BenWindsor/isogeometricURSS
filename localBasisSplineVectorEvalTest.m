U=[0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0];
p=2;
i=3;
hold on;
for j=1:7
    fplot(@(x)(basisSplineEval(U,x,j,p)),[0, 0.999]);
end

point=0.1;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, 1, p));