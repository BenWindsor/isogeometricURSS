% Example from pg119 of 'The NURBS Book'
U=[0,0,0,0,1/4,1/2,3/4,1,1,1,1];
W=[1,1,1,3,1,1,1];
ctrl=[0 1 2 3 5 4.5 6;
      0 2 2.2 -2 -2 1.7 1.8];
p=3;
func=@(x)(nurbCurveEval(U,x,p,ctrl,W));
size=100;
space=linspace(0,0.999,size);
points=zeros(2, size);
for i=1:size
    points(:,i)=func(space(i));
end
scatter(points(1,:), points(2,:));