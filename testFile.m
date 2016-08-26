%General file for testing functions out

U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
p=2;
hold on;
for i=1:7
    fplot(@(x)(basisSplineEval(U,x,i,p)), [0, 0.999])
end
