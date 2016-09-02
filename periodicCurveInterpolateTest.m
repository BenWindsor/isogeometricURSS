%Circle with degree 2 approx
xHandle=@(t)(cos(2*pi*t));
yHandle=@(t)(sin(2*pi*t));

crv=periodicCurveInterpolate(15, 2, xHandle, yHandle);
perbspplot(crv, 100);

figure

%Lissajous curve
xHandle=@(t)(cos(2*pi*3*t));
yHandle=@(t)(sin(2*pi*2*t));

crv=periodicCurveInterpolate(19, 2, xHandle, yHandle);
perbspplot(crv, 300);

% figure
% xHandle=@(t)(t^2);
% crv=periodicCurveInterpolate(20, 2, xHandle);
% perbspplot(crv, 100);

% figure
% xHandle=@(t)(sin(pi*t));
% crv=periodicCurveInterpolate(21, 2, xHandle);
% perbspplot(crv, 100);

%Circle with degree 3 approx
figure
xHandle=@(t)(cos(2*pi*t));
yHandle=@(t)(sin(2*pi*t));

crv=periodicCurveInterpolate(5, 3, xHandle, yHandle);
perbspplot(crv, 100);

