%Computes cg approximations of any Q on time step t_0, t_1
function [soln_mesh_pts, soln_coeffs, soln_times_lagrange, tmesh_lagrange,Y_i_all_imex] = cg_imex_ms(tmesh, Q, y_0, odefe, odegi,other,implicit_table, explicit_table)


N = length(tmesh)-1;

dim = length(y_0);
soln_mesh_pts = zeros(dim,N+1);
soln_coeffs = zeros(N,dim,(Q+1));
soln_mesh_pts(:,1) = y_0;
soln_times_lagrange = zeros(dim,N*Q+1);
soln_times_lagrange(:,1) = y_0;
tmesh_lagrange = ref_mesh(tmesh,Q);

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


imex_coeffs = zeros(dim,nu);

[soln_mesh_pts_imex, Y_i_all_imex, Y_i_ast_all_imex] = imex_ms(tmesh, implicit_table, explicit_table , y_0, odefe, odegi);

for n = 1 : N
    t_0 = tmesh(n);
    t_1 = tmesh(n+1);
    y_init = soln_mesh_pts(:,n);
    other.t_0 = t_0;
    other.t_1 = t_1;
    other.n = n;
    other.N = N;
    
    imex_coeffs(:,:) = Y_i_all_imex(n,:,:);
    [y_coeffs,y_new] = cg_imex_ms_step(t_0, t_1, Q, y_init, odefe,odegi, other,nu,dim,imex_coeffs,w_explicit,w_implicit,c_explicit,c_implicit,A_explicit,A_implicit);
    soln_mesh_pts(:,n+1) = y_new;
    soln_coeffs(n,:,:) = y_coeffs;
    soln_times_lagrange(:, (n-1)*Q+2 : n*Q+1 ) = y_coeffs(:,2:Q+1);
end


end

function [soln,y_1] = cg_imex_ms_step(t_0, t_1, Q, y_0, odefe, odegi, other,nu,dim,imex_coeffs,w_explicit,w_implicit,c_explicit,c_implicit,A_explicit,A_implicit)


%Inputs: y_0: d dimensional vector
%t_0: time at which y_0 is specified
%t_1: final time for the interval. Note that t_0 can be greater than t_1
%Q: what  cg(Q) method to use
%quad_function: provides quad weights and points for integration of (f(y),v)
%odefn: odefn(t,y) is the right hand side of an ode. Since y had
%components, odefn: R^{d+1} -> R

%imex_coeffs = dim x nu sized matrix. imex(:,j) gives the jth stage IMEX variable



%outputs
%soln: coefficients of the basis. (d,(Q+1)) sized matrix. (Q+1) as we also include the first basis, which is known.
%y_1: value of y at time t_1.

% We refer to y_0(i) as the ith component


%%%%%%%% ASSUMING LAGRANGE BASIS FUNCTIONS %%%%%%%%%


dim = length(y_0);


lb = lagrange_basis();
%test_cell_arr has 'Q+1' cells. Each cell is a function of time.
%Note that we know the coefficient for the first basis
test_cell_arr = lb.get_lagrange_basis_funcs(Q+1, t_0, t_1) ;

%test_deriv_cell_arr 'Q+1' cells. Each cell is a function of time.
test_deriv_cell_arr = lb.get_lagrange_basis_funcs_deriv(Q+1, t_0, t_1) ;

%trial_cell_arr has 'Q' cells. Each cell is a function of time.
trial_cell_arr = lb.get_lagrange_basis_funcs(Q, t_0, t_1);





zer = zeros(size(y_0));
%guess = repmat(zer,Q,1);
guess = repmat(y_0,Q,1);

%Now we form the Newton lhs and rhs. Its actually a linear system in this case

res_f = @(curr_soln_vec) Newton_f_res(t_0,t_1, trial_cell_arr, test_cell_arr, test_deriv_cell_arr,Q,dim,odefe, odegi,y_0, curr_soln_vec,other,nu,imex_coeffs,w_explicit,w_implicit,c_explicit,c_implicit,A_explicit,A_implicit);

[soln_vec,ithist,ierr] = nsold(guess, res_f, [1e-14 1e-14], [40 1 0 1]);


%reshape solution
soln = zeros(dim,Q+1);
soln(:,1) = y_0;
for d = 1 : dim
    for q = 1 :Q
        soln(d,q+1) = soln_vec( (q-1)*dim + d);
    end
end

y_1 = soln(:,Q+1);


kkk=1;



end


function [res_vec] = Newton_f_res(t_0,t_1, trial_cell_arr, test_cell_arr, test_deriv_cell_arr,Q,dim,odefe_other,odegi_other,y_0, curr_soln_vec,other,nu,imex_coeffs,w_explicit,w_implicit,c_explicit,c_implicit,A_explicit,A_implicit)

% [ (f_1(y),v_1)
%   (f_2(y),v_1)
%   (f_3(y),v_1)
%     ...
%   (f_{dim}(y),v_1)
%   (f_1(y),v_2)
%     ...

curr_soln = reshape(curr_soln_vec,dim,Q);
y_t_coeff_vec = [y_0 curr_soln];
y_func = @(t) get_y_func_from_basis(t,y_t_coeff_vec,test_cell_arr);
y_func_dot = @(t) get_y_func_from_basis(t,y_t_coeff_vec,test_deriv_cell_arr);
odefe = @(t,y) odefe_other(t,y,other);
odegi = @(t,y) odegi_other(t,y,other);
newton_res_term_tensor = newton_f_residual_term_tensor(t_0,t_1, trial_cell_arr,Q,dim,odefe,odegi, y_func,y_func_dot,nu,imex_coeffs,w_explicit,w_implicit,c_explicit,c_implicit,A_explicit,A_implicit);

res_vec = zeros(Q*dim,1);

for q = 1 : Q
    for d = 1 : dim
        res_vec( (q-1)*dim + d) =  newton_res_term_tensor(d,q);
    end
end

kkk=1;

end

function [integ] = newton_f_residual_term_tensor(t_l,t_r, trial_cell_arr,Q,dim,odefe, odegi, y_func_not_used,y_func_dot,nu,imex_coeffs,w_explicit,w_implicit,c_explicit,c_implicit,A_explicit,A_implicit)

% forms (\dot{y},v)
polydegree = (Q-1)*2;
N_g = ceil(0.5*(polydegree + 1))+1;
integ = zeros(dim,Q);

[x_gauss,w_gauss]=lgwt(N_g,t_l,t_r);

dT = t_r - t_l;

assert( max(abs(c_explicit - c_implicit)) < 1e-14);
x_imex = t_l + c_explicit*(t_r - t_l);

odefe_n = zeros(dim,1);
odegi_n = zeros(dim,1);

for q = 1:Q %all trial functions
    trial_func = trial_cell_arr{q};
    
    
    
    % ydot quadrature
    for n = 1 : N_g %all Gauss points
        for d = 1:dim %all components
            y_comp_deriv = @(t) get_dth_comp(y_func_dot,t,d);
            integ(d,q ) = integ(d,q ) + w_gauss(n)*y_comp_deriv(x_gauss(n))*trial_func(x_gauss(n));
        end
        
    end
    
    
    
    %f(y) quadrature
    for n = 1 : nu
        odefe_n(:) = odefe(x_imex(n),imex_coeffs(:,n));
        odegi_n(:) = odegi(x_imex(n),imex_coeffs(:,n));
        for d = 1:dim %all components
            
            integ(d,q ) = integ(d,q ) - dT*w_explicit(n)*(  odefe_n(d)  *trial_func(x_imex(n)) )...
                - dT*w_implicit(n)*(  odegi_n(d)  *trial_func(x_imex(n)) );
        end
        
    end
   
end



end

function val = get_dth_comp(y,t,d)
arr = y(t);
val = arr(d);

end

function val = get_dim_f_dim_y_comp(y,t,dim_f, dim_y)
mat = y(t);
val = mat(dim_f,dim_y);
end




function tm_refined = ref_mesh(tm_init, factor)
N = length(tm_init)-1;
mesh_pts_new = N*factor + 1;
tm_refined = zeros(mesh_pts_new,1);
tm_refined(1) = tm_init(1);
for n = 1 : N
    dT = tm_init(n+1) - tm_init(n);
    dT_new = dT/factor;
    for r = 1 : factor
        tm_refined( (n-1)*factor + r + 1) = tm_refined( (n-1)*factor + r ) + dT_new;
    end
    tm_refined( n*factor + 1) = tm_init(n+1);
end

kkk=1;
end











