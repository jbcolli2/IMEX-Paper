function [err, I, II, III] = compute_imex_cg_cg_err(soln_coeffs_p,soln_coeffs_d,tmesh_primal,tmesh_dual,primal_Q, dual_Q,p_to_d_dT_ratio,dim,...
    odefe, odegi,implicit_table, explicit_table,Y_i_all_imex)

assert(tmesh_primal(1) == tmesh_dual(1));
assert(tmesh_primal(end) == tmesh_dual(end))


N_dual = length(tmesh_dual) - 1;

err = 0;
I = 0;
II = 0;
III = 0;

y_t_coeff_vec = zeros(dim,primal_Q+1);
phi_t_coeff_vec = zeros(dim,dual_Q+1);


[nup1 dummy] =size(implicit_table);
[dummy_2 dummy_3] =size(explicit_table);
assert(nup1 == dummy_2);
assert(nup1 == dummy_3);
nu = nup1-1;
c_explicit = explicit_table(1:nu,1);
w_explicit  = explicit_table(nup1,2:nup1);
A_explicit = explicit_table(1:nu, 2:nup1);
c_implicit = implicit_table(1:nu,1);
w_implicit  = implicit_table(nup1,2:nup1);
A_implicit = implicit_table(1:nu, 2:nup1);
odefe_n = zeros(dim,1);
odegi_n = zeros(dim,1);
imex_coeffs = zeros(dim,nu);

for n_d = 1 : N_dual
    
    t_l_d = tmesh_dual(n_d);
    t_r_d = tmesh_dual(n_d+1);
    
    n_p = get_primal_time_step_index(N_dual, n_d, p_to_d_dT_ratio,false);
    
    t_l_p = tmesh_primal (n_p);
    t_r_p = tmesh_primal (n_p+1);
    
    lb = lagrange_basis();
    dual_basis_cell_arr = lb.get_lagrange_basis_funcs(dual_Q+1, t_l_d, t_r_d);
    primal_basis_cell_arr = lb.get_lagrange_basis_funcs(primal_Q+1, t_l_p, t_r_p);
    primal_basis_deriv_cell_arr = lb.get_lagrange_basis_funcs_deriv(primal_Q+1, t_l_p, t_r_p);
    
    y_t_coeff_vec_temp = soln_coeffs_p(n_p,:,:);
    y_t_coeff_vec(:,:) = y_t_coeff_vec_temp(1,:,:);
    phi_t_coeff_vec_temp = soln_coeffs_d(n_d,:,:);
    phi_t_coeff_vec(:,:) = phi_t_coeff_vec_temp(1,:,:);
    y = @(t) get_y_func_from_basis(t,y_t_coeff_vec,primal_basis_cell_arr);
    doty = @(t) get_y_func_from_basis(t,y_t_coeff_vec,primal_basis_deriv_cell_arr);
    phi = @(t) get_y_func_from_basis(t,phi_t_coeff_vec,dual_basis_cell_arr);
    
    assert(length(tmesh_dual) == length(tmesh_primal)); %%%Reorder loops if dual_mesh is to be refined
    imex_coeffs(:,:) = Y_i_all_imex(n_p,:,:);
    x_imex = t_l_p + c_explicit*(t_r_p - t_l_p);
    dT = t_r_p - t_l_p;
    
    N_g = 7;
    [x,w]=lgwt(N_g,t_l_d,t_r_d);
    
    odefe_as_t = @(t) odefe(t,y(t));
    odegi_as_t = @(t) odegi(t,y(t));
    
    integ = 0;
    ydotterm  = 0;
    fterm = 0;
    gterm = 0;
    fterm_imex_quad = 0
    gterm_imex_quad = 0;
    for n = 1 : length(x)  %Gauss Quadrature
        t1 =   w(n)* phi(x(n))'*( -doty(x(n)));
        t2 =   w(n)* phi(x(n))'*( odefe_as_t(x(n)));
        t3 =   w(n)* phi(x(n))'*( odegi_as_t(x(n)));
        ydotterm = ydotterm + t1;
        fterm = fterm + t2;
        gterm = gterm + t3;
        integ = integ + t1 + t2 + t3;
    end
    
    for n = 1 : nu %IMEX quadrature
        odefe_n(:) = odefe(x_imex(n),imex_coeffs(:,n));
        odegi_n(:) = odegi(x_imex(n),imex_coeffs(:,n));
        
        
        fterm_imex_quad = fterm_imex_quad + dT*w_explicit(n)*phi(x_imex(n))' *(  odefe_n );
        
        gterm_imex_quad = gterm_imex_quad + dT*w_implicit(n)*phi(x_imex(n))' * (  odegi_n);
        
        
    end
    
    I = I + ydotterm + fterm_imex_quad + gterm_imex_quad;
    II = II  + fterm  - fterm_imex_quad;
    III = III + gterm  - gterm_imex_quad;
    err = err + integ;
    
end
assert( abs(I + II + III - err) < 1e-14);
end
