function rhs = op_f_v_tp_param( space, msh, coeff )
% An analagous version of op_f_v_param.m but taking a function on the
% parameter domain instead of the physical domain. 

 rhs = zeros (space.ndof, 1);


  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);
    
    quadPoints=cell2mat(msh_col.qn);
    rhs = rhs + op_f_v (sp_col, msh_col, coeff (quadPoints));
    
    
  end

end

