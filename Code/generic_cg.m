%Computes cg approximations of any Q on time step t_0, t_1
function [soln_mesh_pts, soln_coeffs, soln_times_lagrange, tmesh_lagrange] = generic_cg(tmesh, Q,quad_function, y_0, odefn,other)





N = length(tmesh)-1;

dim = length(y_0);
soln_mesh_pts = zeros(dim,N+1);
soln_coeffs = zeros(N,dim,(Q+1));
soln_mesh_pts(:,1) = y_0;
soln_times_lagrange = zeros(dim,N*Q+1); 
soln_times_lagrange(:,1) = y_0;
tmesh_lagrange = ref_mesh(tmesh,Q);

for n = 1 : N
    t_0 = tmesh(n);
    t_1 = tmesh(n+1);
    y_init = soln_mesh_pts(:,n);
    other.t_0 = t_0;
    other.t_1 = t_1;
    other.n = n;
    other.N = N;
    [y_coeffs,y_new] = generic_cg_step(t_0, t_1, Q,quad_function, y_init, odefn, other);
    soln_mesh_pts(:,n+1) = y_new;
    soln_coeffs(n,:,:) = y_coeffs;
    soln_times_lagrange(:, (n-1)*Q+2 : n*Q+1 ) = y_coeffs(:,2:Q+1);
end

end

function [soln,y_1] = generic_cg_step(t_0, t_1, Q,quad_function, y_0, odefn, other)


%Inputs: y_0: d dimensional vector
%t_0: time at which y_0 is specified
%t_1: final time for the interval. Note that t_0 can be greater than t_1
%Q: what  cg(Q) method to use
%quad_function: provides quad weights and points for integration of (f(y),v)
%odefn: odefn(t,y) is the right hand side of an ode. Since y had
%components, odefn: R^{d+1} -> R
%odefn_d_y(t,y): Jacobian of odefn. R^{d+1} -> R^d x R^d

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

%Now we form the Newton lhs and rhs

res_f = @(curr_soln_vec) Newton_f_res(t_0,t_1, trial_cell_arr, test_cell_arr, test_deriv_cell_arr,Q,dim,odefn,quad_function,y_0, curr_soln_vec,other);

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


function [res_vec] = Newton_f_res(t_0,t_1, trial_cell_arr, test_cell_arr, test_deriv_cell_arr,Q,dim,odefn_w_other,quad_function,y_0, curr_soln_vec,other)

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
odefn = @(t,y) odefn_w_other(t,y,other);
newton_res_term_tensor = newton_f_residual_term_tensor(t_0,t_1, trial_cell_arr,Q,dim,odefn, y_func,y_func_dot,quad_function);

res_vec = zeros(Q*dim,1);

for q = 1 : Q
    for d = 1 : dim
        res_vec( (q-1)*dim + d) =  newton_res_term_tensor(d,q);
    end
end

kkk=1;

end

function [integ] = newton_f_residual_term_tensor(t_l,t_r, trial_cell_arr,Q,dim,odefn, y_func,y_func_dot,quad_function)

% forms (\dot{y},v)
polydegree = (Q-1)*2;
N_g = ceil(0.5*(polydegree + 1))+1;
integ = zeros(dim,Q);

[x_gauss,w_gauss]=lgwt(N_g,t_l,t_r);

[x_odefunc,w_odefunc]=quad_function(t_l,t_r);

for d = 1:dim %all components
    %y_comp = @(t) get_dth_comp(y_func,t,d);
    y_comp_deriv = @(t) get_dth_comp(y_func_dot,t,d);
    odefn_as_t = @(t) odefn(t,y_func(t));
    odefn_comp = @(t) get_dth_comp(odefn_as_t,t,d);
    for q = 1:Q %all trial functions
        trial_func = trial_cell_arr{q};
        
        % ydot quadrature
        for n = 1 : N_g %all Gauss points
            integ(d,q ) = integ(d,q ) + w_gauss(n)*y_comp_deriv(x_gauss(n))*trial_func(x_gauss(n));
        end
        
        %f(y) quadrature
        for n = 1 : length(x_odefunc)
            integ(d,q ) = integ(d,q ) - w_odefunc(n)*odefn_comp(x_odefunc(n))*trial_func(x_odefunc(n));
        end
        
    end
end



end



function val = get_dth_comp(y,t,d)
arr = y(t);
val = arr(d);

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











