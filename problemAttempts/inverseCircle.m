function point = inverseCircle( x, y)
%Take points x,y on circle and return the point on the unit inverval that
%would have been mapped to it by the map f(x) = (cos(2Pix) , sin(2pix))


if numel(x)==numel(y)
    point=zeros(numel(x),1);
else
    error('Need same number of x and y vals');
end

for j=1:numel(point)
    %Normalise/project to circle if dealing with approximate circle
    size=sqrt(x(j)*x(j) + y(j)*y(j));
    xNew=x(j)/size;
    yNew=y(j)/size;
    if yNew>=0
        point(j) = (1/(2*pi))*acos(xNew);
    else
        %point(j) = (1/(2*pi))*(acos(xNew)+pi);
        point(j)=1-(1/(2*pi))*acos(xNew);
    end
end


end

