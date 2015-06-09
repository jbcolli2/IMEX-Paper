classdef imex_stuff
    
    methods
        
        function [implicit_table, explicit_table,primal_Q] =  get_mid_point_tables(obj)
            implicit_table = [0,0 ,0; 0.5 0 0.5; nan 0 1];
            explicit_table = [0,0 ,0; 0.5 0.5 0; nan 0 1];
            primal_Q = 1;
        end

        
        
    end
    
    
end