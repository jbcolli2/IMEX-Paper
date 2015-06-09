classdef ode_problem
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        
    end
    
    methods
        function o = ode_problem()
            
        end
        
     
  
        function val = get_u_true(obj,t)
            'Computing U True'
            %options=odeset('RelTol',1e-12);
            options=odeset('RelTol',1e-12,'AbsTol',1e-12);
            Xo = obj.get_init();
            %timespan
            tspan = [0,t];
            %call the solver
            TestFunction = @(t,u) obj.ode_fun(t,u);
            [t,X] = ode23s(TestFunction,tspan,Xo,options);
            nn = length(t);
            val = X(nn,:)';
            
        end
        
        function [t,X] = get_u_true_full(obj,t)
            'true val'
            options=odeset('RelTol',3e-14,'AbsTol',3e-14);           
            Xo = obj.get_init();
            %timespan
            tspan = [0,t];
            %call the solver            
            TestFunction = @(t,u) obj.ode_fun(t,u);
            [t,X] = ode23s(TestFunction,tspan,Xo,options);
           
        end
        
        function val =  neg_f_prime_transpose_times_phi(obj,t,phi,other) % neg so that we only have \dot{\phi} on lhs
            
            %other has information about the primal solution
            y_coeff_tensor =  other.y_coeff_tensor;
            p_to_d_dT_ratio = other.p_to_d_dT_ratio;
            t_0 = other.t_0 ;
            t_1 = other.t_1;
            n_dual = other.n;
            N_dual = other.N;
            t_mesh_primal = other.tmesh_primal;
            primal_Q = other.primal_Q;
            
            assert(t_1 < t_0); % we are solving duals
            
            n_primal = get_primal_time_step_index(N_dual, n_dual, p_to_d_dT_ratio,true);
            
            t_l = t_mesh_primal (n_primal);
            t_r = t_mesh_primal (n_primal+1);
            lb = lagrange_basis();
            primal_basis_cell_arr = lb.get_lagrange_basis_funcs(primal_Q+1, t_l, t_r);
            
            assert(t_r >= t_0);
            assert(t_l <= t_1);
            assert( t <= t_r);
            assert(t >= t_l);
            
            
            y_t_coeff_vec_temp = y_coeff_tensor(n_primal,:,:);
            assert(ndims(y_t_coeff_vec_temp) == 3);
            [dummy, dim,num_basis] = size(y_t_coeff_vec_temp);
            
            y_t_coeff_vec = zeros(dim,num_basis);
            y_t_coeff_vec(:,:) = y_t_coeff_vec_temp(1,:,:);
            u = get_y_func_from_basis(t,y_t_coeff_vec,primal_basis_cell_arr);
            
            
            val =  (-1.*obj.d_ode_fun_u(t,u))'*phi; % -1 so that we only have \dot{\phi} on lhs
            
        end
          
        function [T,dT,dT_d_p_ratio,dual_order_increment]  = get_suggested_parms(obj)
            T = [];
            dT = [];
            dT_d_p_ratio = [];
            dual_order_increment = 1;
        end
          
    end
    
     methods (Abstract)
        get_name(obj);
        ode_fun(obj,t,u,other);
        d_ode_fun_u(obj,t,u);
        get_dim(obj);
        get_init(obj);
        get_psi_N(obj);
        
     end
end