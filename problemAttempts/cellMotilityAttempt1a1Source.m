function vals = cellMotilityAttempt1a1Source( preva1, preva2, preva3, gamma, r1, s, s1, s3, b1, x )

vals=zeros(1,numel(x));

for i=1:numel(x)
    vals(i)=gamma*((((r1+s)*(perbspeval(preva1,x(i))*perbspeval(preva1,x(i))/preva2 + b1))/((s3+perbspeval(preva3, x(i)))*(1+s1*perbspeval(preva1,x(i))*perbspeval(preva1,x(i)))))-r1*perbspeval(preva1,x(i)));
end

end

