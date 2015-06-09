classdef sim_parms
    
    properties
        
        T = 5.0;
        dT_d_p_ratio = 1;
        dT = 0.01;
        dual_order_increment = 0;
    end
    
    methods
        function obj = set_suggested_parms(obj,T_sug,dT_sug,dT_d_p_ratio_sug,dual_order_increment_sug)
            if ~isempty(T_sug)
                obj.T  = T_sug;
            end
            if ~isempty(dT_sug)
                obj.dT = dT_sug;
            end
            if ~isempty(dT_d_p_ratio_sug)
                obj.dT_d_p_ratio = dT_d_p_ratio_sug;
            end
            if ~isempty(dual_order_increment_sug)
                obj.dual_order_increment = dual_order_increment_sug;
            end
            
        end
    end
end