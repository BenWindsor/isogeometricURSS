%1D test
U = [0 0.1 0.2 0.3 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
W = [1 0.2 1 1 0.2 1 1 0.1];
p = 2;
ctrl=[0 1 2.5 2.6 2.7 2.8 6 2];
func=@(x)(nurbCurveEval(U, x, p, ctrl, W));
fplot(func, [0 1])

%2D test
% U = [0,0,0,1,1,1];
% W=[1, 1/sqrt(2), 1];
% p = 1;
% ctrl=[0 -1 -1;
%       1 1 0];
% func = @(x)(nurbCurveEval(U,x,p,ctrl,W));
% size=50;
% space=linspace(0,1,size);
% points=zeros(2, size);
% for i=1:size
%     points(:,i)=func(space(i));
% end
% scatter(points(1,:), points(2,:));


