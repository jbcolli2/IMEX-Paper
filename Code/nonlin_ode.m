classdef nonlin_ode < imex_ode_problem
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties        
        lambda = -2;
        
    end
    
    methods
        function o = nonlin_ode()
             o = o@imex_ode_problem();
        end
    end
    
     methods 
         function s = get_name(obj)
             s = 'nonlin_ode';
         end
         
        
         function val = f(obj,t,u,other)
            val =   .1*(u*u-16);
         end
         
         function val = g(obj,t,u,other)
            val =  exp(.1*u).*sin(u) ;
         end

        function val = df_u(obj,t,u)
            val =  .1*(2*u);
            
        end
         
        function val = dg_u(obj,t,u)
            val = .1*exp(.1*u).*sin(u) + exp(.1*u).*cos(u);
            
         end
        
        function val = get_init(obj)
            val = 10;
        end
        
        function val = get_psi_N(obj)
            val = 1;
        end
        
        function dim =  get_dim(obj)
            dim = 1;
        end
        
        function [T,dT,dT_d_p_ratio,dual_order_increment]  = get_suggested_parms(obj)
            [T,dT,dT_d_p_ratio,dual_order_increment] = get_suggested_parms@imex_ode_problem(obj);
            T = 0.8;
            dT = 0.01;
            dT_d_p_ratio = 1;
        end
        
       
     end
end