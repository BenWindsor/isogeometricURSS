% U=[0 0 0 1 2 3 3 3];
% ctrl=[0 0.25 0.75 1 1;
%        0 1 0.6 0.2 0];
% p=2;
% func=@(x)(periodicSplineCurveEval(U,x,p,ctrl));
% points=func(linspace(0,2.999,100));
% scatter(points(1,:), points(2,:));

%Plot a 1D curve
U=[0 0 0 0.2 0.4 0.6 0.8 1 1 1];
ctrl=[1 2 3 4 -3];
p=2;
fplot(@(x)(periodicSplineCurveEval(U,x,p,ctrl)), [0, 0.999]);