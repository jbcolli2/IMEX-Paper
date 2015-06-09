classdef lagrange_basis
    
    methods
        function basis_derivs = get_lagrange_basis_funcs_deriv(self,num_basis, t_0, t_1)
            basis_derivs = cell(num_basis,1);
            
            if num_basis == 1
                basis_derivs{1} = @(t) zeros(size(t));
            else
                dT = (t_1 - t_0)/(num_basis-1);
                pts = t_0:dT:t_1;
                if length(pts) ~= num_basis
                    kkk=1;
                end
                assert(length(pts) == num_basis);
                pts(end) = t_1;
                for i = 1 : num_basis
                    basis_derivs{i}  = @(t)self.get_ith_lagrange_basis_deriv(t,pts,i);
                end
            end
        end
        
        function val = get_ith_lagrange_basis_deriv(self,t,pts,i)
            num_pts = length(pts);
            den = ones(size(t));
            
            for p = 1 : num_pts
                if p ~= i
                    den = den* (pts(i) - pts(p))  ;
                end
            end
            
            num = zeros(size(t));
            for k = 1 : num_pts
                if k ~= i
                    num_summand_term = ones(size(t));
                    for p = 1 : 1:num_pts
                        if p ~= i && p ~= k
                            num_summand_term = num_summand_term.*(t - pts(p));
                        end
                    end
                    num = num + num_summand_term;
                end
            end
            
            val = num./den;
        end
        
        
        function basis = get_lagrange_basis_funcs(self,num_basis, t_0, t_1)
            basis = cell(num_basis,1);
            
            if num_basis == 1
                basis{1} = @(t) ones(size(t));
            else
                dT = (t_1 - t_0)/(num_basis-1);
                pts = t_0:dT:t_1;
                assert(length(pts) == num_basis);
                pts(end) = t_1;
                for i = 1 : num_basis
                    basis{i}  = @(t)self.get_ith_lagrange_basis(t,pts,i);
                end
            end
            
            
            
            
        end
        
        function val = get_ith_lagrange_basis(self,t,pts,i)
            num_pts = length(pts);
            val = ones(size(t));
            for p = 1 : num_pts
                if p ~= i
                    val = val.* ( (t - pts(p)) ./ ((pts(i) - pts(p)))  );
                end
            end
        end
        
    end
    
end