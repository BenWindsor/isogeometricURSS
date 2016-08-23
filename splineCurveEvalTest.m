% 1D test
% U = [0 0.1 0.2 0.3 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
% p = 2;
% ctrl=[0 1 2.5 2.6 2.7 2.8 6 7];
% func=@(x)(splineCurveEval(U, x, p, ctrl));
% fplot(func, [0 1])

%2D test
U = [0 0.1 0.2 0.3 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
p = 2;
ctrl=[0 1 2.5 2.6 2.7 2.8 6 7;
      0 1 2.2 2.4 2.5 2.2 1 0];
func = @(x)(splineCurveEval(U,x,p,ctrl));
points=func(linspace(0,1,500));
scatter(points(1,:), points(2,:));

