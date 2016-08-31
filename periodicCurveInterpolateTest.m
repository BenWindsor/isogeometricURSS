%Circle
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
% xHandle=@(t)(cos(pi*t));
% crv=periodicCurveInterpolate(20, 2, xHandle);
% perbspplot(crv, 100);