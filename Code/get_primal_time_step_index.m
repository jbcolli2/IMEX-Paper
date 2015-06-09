function n_primal = get_primal_time_step_index(N_dual, n_dual, p_to_d_dT_ratio,is_time_reverse)
%assumes that the mesh is uniformly refined

if is_time_reverse
    n_dual = N_dual - n_dual + 1;
end
n_primal = n_dual/p_to_d_dT_ratio;
if ( abs(n_primal - round(n_primal)) < 1e-10)
    n_primal = round(n_primal);
else
    n_primal = ceil(n_primal);
end