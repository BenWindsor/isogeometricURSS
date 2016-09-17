function vals = cellMotilityAttempt1AdvectionField( prevCrv, newCrv, delta, x)
% Return the advection field of the curve at points x

vals=zeros(2, numel(x));

for i=1:numel(x)
    %vals(i) = periodicSplineCurveDerivEval(prevCrv, x(i));
    vals(:,i) = (perbspeval(newCrv, x(i))-perbspeval(prevCrv, x(i)))/delta;
end

end

