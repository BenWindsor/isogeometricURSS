%plot the local functions on [0, 0.2)
U=[0 0 0 0 0.2 0.4 0.6 0.8 1.0 1.0 1.0 1.0];
p=3;
i=1;
size=100;
space=linspace(U(4), U(4+1), size);
points=zeros(size, p+1);
for j=1:numel(space)
    points(j,:)=localBasisSplineVectorEval(U,space(j),4,p);
end
hold on;
for j=1:(p+1)
    scatter(space, points(:,j));
end

%Convert the points with the operator and re plot
newPoints=points;
operator=localPeriodicOperatorDeg3(U,1);
for j=1:numel(space)
    point=transpose(newPoints(j,:));
    newPoints(j,:)=operator*point;
end
for j=1:(p+1)
    scatter(space, newPoints(:,j), 'filled');
end

%Print next section
space=linspace(U(5), U(5+1), size);
points=zeros(size, p+1);
for j=1:numel(space)
    points(j,:)=localBasisSplineVectorEval(U,space(j),5,p);
end
for j=1:(p+1)
    scatter(space, points(:,j));
end