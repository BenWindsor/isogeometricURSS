function pernurb = perbspmak( coefs, knots )
% The analogue of the nurbs toolbox nrbmak.m function producing a periodic
% nurb structure
pernurb = struct ('form', 'Per-BSpline', 'dim', 2, 'number', [], 'coefs', [], ...
                'knots', [], 'order', []);            
pernurb.form   = 'Per-BSpline';
pernurb.dim    = 2;
np = size(coefs);
dim = np(1);
pernurb.number=np(2); %number of control points
pernurb.coefs=coefs;
pernurb.order=size(knots,2)-np(2); %one more than degree
knots=sort(knots);
pernurb.knots=knots;

end

