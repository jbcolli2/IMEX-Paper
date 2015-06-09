classdef imex_scalar_ode < imex_ode_problem
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        lambda = -2;
        
    end
    
    methods
        function o = imex_scalar_ode()
             o = o@imex_ode_problem();
        end
    end
    
     methods 
         function s = get_name(obj)
             s = 'Scalar ODE';
         end
         
        function val = f(obj,t,u,other)
            val = 0.5*obj.lambda * u;
        end
        function val = g(obj,t,u,other)
            val = 0.5*obj.lambda * u;
        end
        function val = df_u(obj,t,u)
            %assert(0);
            val = 0.5*obj.lambda *ones(size(t));
        end
        function val = dg_u(obj,t,u)
            
            val = 0.5*obj.lambda *ones(size(t));
        end
        function val = get_init(obj)
            val = 10;
        end
        function val = get_u_true(obj,t)
            val = obj.get_init()*exp(obj.lambda*t);
        end
        function val = get_psi_N(obj)
            val = 1;
        end
        
        function dim =  get_dim(obj)
            dim = 1;
        end
        
        function [T,dT,dT_d_p_ratio,dual_order_increment]  = get_suggested_parms(obj)
            T = 0.8;
            dT = 0.1;
            dT_d_p_ratio = 1;
            dual_order_increment = 1;
        end
        
       function val =  neg_f_prime(obj,t,other) % neg so that we only have \dot{\phi} on lhs
            dh_u = @(t,u) obj.df_u(t,u) + obj.dg_u(t,u);
            val =  -1.*obj.lambda *ones(size(t)); % -1 so that we only have \dot{\phi} on lhs
        end
     end
end