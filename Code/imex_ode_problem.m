classdef imex_ode_problem < ode_problem
    
    properties
        
        
    end
    
    methods
        function o = imex_ode_problem()
            
        end
  
        function val = ode_fun(obj,t,u,other)
           val = obj.f(t,u) + obj.g(t,u);
        end
        
        function val =  d_ode_fun_u(obj,t,u)
            val = obj.df_u(t,u) + obj.dg_u(t,u);
        end
    end
    
     methods (Abstract)
        
        
        f(obj,t,u);
        g(obj,t,u);
        df_u(obj,t,u);
        dg_u(obj,t,u);
        
     end
end