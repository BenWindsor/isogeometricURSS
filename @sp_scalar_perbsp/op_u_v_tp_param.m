function varargout = op_u_v_tp_param (space1, space2, msh, coeff)
% Analogue of op_u_v_tp.m but with coeff defined on the parametric domain
% rather than on the mesh/curve, so no pullback required. 

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col);
    sp2_col = sp_evaluate_col (space2, msh_col);

    if (nargin == 4)
      quadPoints=cell2mat(msh_col.qn);
      coeffs = coeff (quadPoints);
    else
      coeffs = ones (msh_col.nqn, msh_col.nel);
    end

    A = A + op_u_v (sp1_col, sp2_col, msh_col, coeffs);
  end

  % As there are less periodic splines the last two rows and cols are zero
  % so cut them off
  sp_univ1=space1.sp_univ;
  sp_univ2=space2.sp_univ;
  space1Size=size(sp_univ1.connectivity);
  space2Size=size(sp_univ2.connectivity);
  elemNum1=space1Size(2);
  elemNum2=space2Size(2);
  A=A(1:elemNum1, 1:elemNum2);
  
  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
