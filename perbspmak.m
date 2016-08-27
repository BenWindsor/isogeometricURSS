function pernurb = perbspmak( coefs, knots )
% The analogue of the nurbs toolbox nrbmak.m function producing a periodic
% nurb structure
pernurb = struct ('form', 'Per-BSpline', 'dim', 2, 'number', [], 'coefs', [], ...
                'knots', [], 'order', []);            
pernurb.form   = 'Per-BSpline';
np = size(coefs);
dim = np(1);
pernurb.dim    = dim;
pernurb.number=np(2); %number of control points
pernurb.coefs=coefs;
%pernurb.order=size(knots,2)-np(2); %one more than degree
pernurb.order=3;
knots=sort(knots);
pernurb.knots=knots;

end

