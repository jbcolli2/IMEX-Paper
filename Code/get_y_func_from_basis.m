function y_func = get_y_func_from_basis(t,y_t_coeff_vec,basis_cell_arr)

if ndims(y_t_coeff_vec) ~= 2
    kkk=1;
end
assert(ndims(y_t_coeff_vec) == 2);

[dim,num_basis] = size(y_t_coeff_vec);
y_func = zeros(dim,1);
for d = 1 : dim
    coeff_d = y_t_coeff_vec(d,:);
    for q = 1 : num_basis
        basis = basis_cell_arr{q};
        y_func(d) = y_func(d) + coeff_d(q)*basis(t);
        
    end
end


end
