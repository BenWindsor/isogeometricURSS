%Plot some basis functions
U = [0,0,0,0,0.33,0.67,1,1,1,1];
p=3;
W=[1 1 1 1 1 1];
hold on;
for i=1:5
    fplot(@(x)(nurbBasisEval(U,x,i,p,W)), [0 0.9999]);
end
%Plot the returned points to see if they fall where they should on the
%functions
point=0.5;
scatter(point*ones(p+1,1), localNurbVectorEval(U, point, 5, p, W));
point=0.3;
scatter(point*ones(p+1,1), localNurbVectorEval(U, point, 5, p, W));
