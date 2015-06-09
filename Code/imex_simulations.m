function imex_simulations()

% Setup. ODE Problem and parameters
my_ode = nonlin_ode();
%my_ode = imex_scalar_ode();

sparms = sim_parms();
imex_tables = imex_stuff();
[T_sug,dT_sug,dT_d_p_ratio_sug,dual_order_increment]  = my_ode.get_suggested_parms();
sparms = sparms.set_suggested_parms(T_sug,dT_sug,dT_d_p_ratio_sug,dual_order_increment);

y_init=my_ode.get_init();
odefe =@my_ode.f;
odegi =@my_ode.g;

[implicit_table, explicit_table,primal_Q] =  imex_tables.get_mid_point_tables();
other = [];


%%%%%%FORWARD (or PRIMAL)%%%%%%
tmesh = [0:sparms.dT:sparms.T];
[soln_mesh_pts_p, soln_coeffs_p, soln_times_lagrange_p, tmesh_lagrange_p,Y_i_all_imex] = cg_imex_ms(tmesh, primal_Q, y_init, odefe, odegi,other,implicit_table, explicit_table);



%%%%%%ADJOINT%%%%%%
other = other_parms();
other.y_coeff_tensor = soln_coeffs_p;

% The adjoint mesh could be a uniform refinement by a factor of dT_d_p_ratio
% of the primal mesh
other.p_to_d_dT_ratio = sparms.dT_d_p_ratio;

other.tmesh_primal = tmesh;
other.primal_Q = primal_Q;
myquads = my_quads();
quad_function = @(t_l,t_r)myquads.gauss_quad(7,t_l,t_r);

% The adjoint basis is higher than the forward basis by dual_order_increment
dual_Q = primal_Q + sparms.dual_order_increment;

phi_0 = my_ode.get_psi_N();
odefn = @my_ode.neg_f_prime_transpose_times_phi;
tmesh_dual = tmesh(end:-1:1);
[soln_mesh_pts_d_r, soln_coeffs_d_r, soln_times_lagrange_d_r, tmesh_lagrange_d_r] = generic_cg(tmesh_dual, dual_Q,quad_function, phi_0, odefn, other);


%Reverse dual coefficients so that they are also "Forward" in time
soln_coeffs_d = zeros(size(soln_coeffs_d_r));
N_dual = length(tmesh_dual) - 1;
for n = 1 : N_dual
   k = N_dual - n + 1;
   soln_coeffs_d(k,:,1:dual_Q+1) = soln_coeffs_d_r(n,:,dual_Q+1:-1:1);
end


%Compute Error
odefn=@my_ode.ode_fun;
dim = length(y_init);
tmesh_dual = tmesh_dual(end:-1:1);


%err_comput_old = compute_cg_cg_err(soln_coeffs_p,soln_coeffs_d,tmesh,tmesh_dual,primal_Q, dual_Q,sparms.dT_d_p_ratio, odefn,dim)

format long E;

[err_comput, I, II, III] = compute_imex_cg_cg_err(soln_coeffs_p,soln_coeffs_d,tmesh,tmesh_dual,primal_Q, dual_Q,sparms.dT_d_p_ratio,dim,odefe, odegi,...
    implicit_table, explicit_table,Y_i_all_imex)

err_actual =   my_ode.get_psi_N()'*(my_ode.get_u_true(sparms.T) - soln_mesh_pts_p(:,end))





end