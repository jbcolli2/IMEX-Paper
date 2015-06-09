classdef my_quads
    
    methods
        function [x,w] = crank_nicolson(obj,t_l,t_r) 
            x = [t_l, t_r];
            dT = t_r - t_l;
            w = [0.5*dT, 0.5*dT];
            
        end
        function [x,w] = gauss_quad(obj,N,t_l,t_r) 
           [x,w]=lgwt(N,t_l,t_r); 
        end
        
    end
    
    
end