function sp = sp_perbsp_1d_param( knots, degree, nodes, varargin )
% An analagous version of sp_bspline_1d_param.m but for periodic basis
% splines

%IMPLEMENT: handle hessian, so set up second derivatives.

%Check arguments in to see if gradient specified true or false. True by
%default
gradient = true; 
if (~isempty (varargin))
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    else 
      error ('sp_bspline_1d_param: unknown option %s', varargin {ii});
    end
  end
end

%Assign number of derivatives, nders, depending on options given. 
if(gradient)
  nders=1;
else
  nders=0;
end

% set up properties
mknots = length (knots)-1; %Chainging to start from zero, MAYBE DONT DO THIS SO IT WILL WORK WITH MY FIND SPAN...
p      = degree;
mcp    = -p - 1 + mknots; %ALTER: check this is correct for control points number in periodic splines
ndof   = mcp + 1;

%Set number of elements and quadrature nodes. 
nel = size (nodes, 2);
nqn = size (nodes, 1);

%Set up connectivity array
nsh = zeros (1, nel);
connectivity = zeros (p+1, nel);
for iel=1:nel
  %s = findspan(mcp, p, nodes(:, iel)', knots);
  s = findSpan(knots,nodes(:,iel)'); 
  %c = numbasisfun(s, nodes(:, iel)', p, knots); %WRITE numperbasisfun.m
  c = numperbsp(knots, nodes(:,iel)', p);
  %c = unique(c(:))+1; % QUESTION: do I remove the +1 here?? as we are now counting from 0 anyway?
  c = unique(c(:));
  connectivity(1:numel(c), iel) = c;
  nsh(iel) = nnz (connectivity(:,iel));
end

nsh_max = max (nsh);

% s = findspan(mcp, p, nodes, knots);
s = findSpan(knots, nodes);
%tders = basisfunder (s, p, nodes, knots, nders); %ALTER to use periodic derivative
tders = perbspfunder(knots, nodes, p, nders);
%nbf   = numbasisfun (s(:)', nodes(:)', p, knots); %ALTER to use numperbasisfun.m
nbf = numperbsp(knots, nodes(:)', p);
%nbf   = reshape (nbf+1, size(s,1), size(s,2), p+1);
nbf = reshape(nbf, size(s,q), size(s,2), p+1);

ders = zeros (numel(nodes), nders+1, nsh_max);
for inqn = 1:numel(nodes)
  [ir,iel] = ind2sub (size(nodes),inqn);
  ind = find (connectivity(:,iel) == nbf(ir,iel,1)); 
  ders(inqn,:,ind:ind+p) = tders(inqn,:,:);
end

supp = cell (ndof, 1);
for ii = 1:ndof
  [dummy, supp{ii}] = find (connectivity == ii);
end

shape_functions = reshape (ders(:, 1, :), nqn, nel, []);
shape_functions = permute (shape_functions, [1, 3, 2]);

sp = struct ('nsh_max', nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
             'connectivity', connectivity, ...
             'shape_functions', shape_functions, ...
             'ncomp', 1, 'degree', degree, 'knots', knots);
sp.supp = supp;

if (gradient)
  shape_function_gradients = reshape (ders(:, 2, :), nqn, nel, []);
  shape_function_gradients = permute (shape_function_gradients, [1, 3, 2]);
  sp.('shape_function_gradients') = shape_function_gradients;
end


end

