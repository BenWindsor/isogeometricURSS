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
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    else 
      error ('sp_bspline_1d_param: unknown option %s', varargin {ii});
    end
  end
end

%Assign number of derivatives, nders, depending on options given. 
if (hessian)
  nders=2;
elseif(gradient)
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
    %   c = numperbsp(knots, nodes(:,iel)', p);
    %c = unique(c(:))+1; % QUESTION: do I remove the +1 here?? as we are now counting from 0 anyway?
    %   c = unique(c(:));
    
    if p==2
        %Altered to make correct connectivity array by hand
        if iel==1
            connectivity(1:3, iel) = [nel-1 nel 1];
        else
            if iel==2
                connectivity(1:3, iel) = [nel 1 2];
            else
                connectivity(1:3, iel) = [iel-2 iel-1 iel];
            end
        end
    elseif p==3
        if iel==1
            connectivity(1:4, iel) = [nel-2 nel-1 nel 1];
        else
            if iel==2
                connectivity(1:4, iel) = [nel-1 nel 1 2];
            elseif iel==3
                connectivity(1:4, iel) = [nel 1 2 3];
            else
                connectivity(1:4, iel) = [iel-3 iel-2 iel-1 iel];
            end
        end
    end
    %  connectivity(1:numel(c), iel) = c;
    nsh(iel) = nnz (connectivity(:,iel));
end

nsh_max = max (nsh);

% s = findspan(mcp, p, nodes, knots);
s = findSpan(knots, nodes);
%tders = basisfunder (s, p, nodes, knots, nders); 
tders = perbspfunder(knots, nodes, p, nders);
%nbf   = numbasisfun (s(:)', nodes(:)', p, knots); 
nbf = numperbsp(knots, nodes(:)', p);
%nbf   = reshape (nbf+1, size(s,1), size(s,2), p+1);
nbf = reshape(nbf, size(s,1), size(s,2), p+1);


ders = zeros (numel(nodes), nders+1, nsh_max);
%DEBUG  --------------------------------------------------------
% fprintf('sp_perbsp_1d_param.m debug line 65\n');
% fprintf('\n');
% fprintf('size(ders) before: ');
% fprintf(mat2str(size(ders)));
% fprintf('\n');
%END DEBUG------------------------------------------------------

for inqn = 1:numel(nodes)
  [ir,iel] = ind2sub (size(nodes),inqn);
  %TEST FIX for resize changing nunber of elements bug
  %ind = find (connectivity(:,iel) == nbf(ir,iel,1)); 
  ind=1;
  %END TEST FIX
  ders(inqn,:,ind:ind+p) = tders(inqn,:,:);
    %DEBUG  --------------------------------------------------------
%     fprintf('\n');
%     fprintf('size(ders) during loop: ');
%     fprintf(mat2str(size(ders)));
%     fprintf('\n');
%     fprintf('ind: ');
%     fprintf(mat2str(ind));
%     fprintf('\n');
    %END DEBUG------------------------------------------------------
end

% DEBUG ---------------------------------------------------------
% fprintf('sp_perbsp_1d_param.m debug line 84\n');
% fprintf('tders: ');
% fprintf(mat2str(tders(:,:,1)));
% % fprintf('nbf: ');
% % fprintf(mat2str(nbf));
% fprintf('numel(nodes): ');
% fprintf(num2str(numel(nodes)));
% fprintf('\n');
% fprintf('nders+1: ');
% fprintf(num2str(nders+1));
% fprintf('\n');
% fprintf('nsh_max: ');
% fprintf(num2str(nsh_max));
% fprintf('\n');
% fprintf('ders(:,:,1): ');
% fprintf(mat2str(ders(:,:,1)));
% fprintf('\n');
% fprintf('size(ders) after: ');
% fprintf(mat2str(size(ders)));
% fprintf('\n');
% fprintf('nsh: ');
% fprintf(mat2str(nsh));
% fprintf('\n');
% END DEBUG --------------------------------------------------------

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

%DEBUG-------------------------------------------------------
% fprintf('shape_functions(:,:,1): ');
% fprintf(mat2str(shape_functions(:,:,1)));
% fprintf('\n');
%END DEBUG -----------------------------------------------------

if (gradient)
  shape_function_gradients = reshape (ders(:, 2, :), nqn, nel, []);
  shape_function_gradients = permute (shape_function_gradients, [1, 3, 2]);
  sp.('shape_function_gradients') = shape_function_gradients;
end

if (hessian)
  shape_function_hessians = reshape (ders(:, 3, :), nqn, nel, []);
  shape_function_hessians = permute (shape_function_hessians, [1, 3, 2]);
  sp.('shape_function_hessians') = shape_function_hessians;
end

end

