function vals = cellMotilityAttempt1a3Source( preva1, preva3, gamma, r3, b3, x )

vals=zeros(1, numel(x));

for i=1:numel(x)
    vals(i)=gamma*(b3*perbspeval(preva1, x(i)) - r3*perbspeval(preva3, x(i)));
end

end

