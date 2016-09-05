
% Plot some curves and points for degree 1
U=[0 0 0.2 0.4 0.6 0.8 1.0 1.0];
p=1;
i=3;
hold on;
for j=1:6
    fplot(@(x)(basisSplineEval(U,x,j,p)),[0, 0.999]);
end

%Plot the points vector at a specific value and make sure they land on the
%curves
point=0.09;
elem=1;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

point=0.45;
elem=3;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

point=0.95;
elem=5;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));


figure
% Same for degree 2
U=[0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0];
p=2;
i=3;
hold on;
for j=1:7
    fplot(@(x)(basisSplineEval(U,x,j,p)),[0, 0.999]);
end

%Plot the points vector at a specific value and make sure they land on the
%curves
point=0.09;
elem=1;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

point=0.45;
elem=3;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

point=0.95;
elem=5;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

% Same for degree 3
figure
%Plot some basis curves
U=[0 0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0 1.0];
p=3;
i=3;
hold on;
for j=1:8
    fplot(@(x)(basisSplineEval(U,x,j,p)),[0, 0.999]);
end

%Plot the points vector at a specific value and make sure they land on the
%curves
point=0.09;
elem=1;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

point=0.45;
elem=3;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));

point=0.95;
elem=5;
scatter(point*ones(p+1,1), localBasisSplineVectorEval(U, point, elem, p));