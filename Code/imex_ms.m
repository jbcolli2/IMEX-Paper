function [soln_mesh_pts, Y_i_all, Y_i_ast_all] = imex_ms(tmesh, implicit_table, explicit_table , y_init, odef, odeg)

[nup1 dummy] =size(implicit_table);
[dummy_2 dummy_3] =size(explicit_table);
assert(nup1 == dummy_2);
assert(nup1 == dummy_3);
nu = nup1-1;
c_explicit = explicit_table(1:nu,1);
w_explicit  = explicit_table(nup1,2:nup1);
A_explicit = explicit_table(1:nu, 2:nup1);
c_implicit = implicit_table(1:nu,1);
w_implicit  = implicit_table(nup1,2:nup1);
A_implicit = implicit_table(1:nu, 2:nup1);

dim = length(y_init);
N = length(tmesh)-1;

Y_i_ast = zeros(dim,nu);
Y_i = zeros(dim,nu);
soln_mesh_pts = zeros(dim,N+1);
soln_mesh_pts(:,1) = y_init;
Y_i_all = zeros(N,dim,nu);
Y_i_ast_all = zeros(N,dim,nu);

for n = 1 : N
    t_0 = tmesh(n);
    t_1 = tmesh(n+1);
    dT = t_1 - t_0;
    y_0 = soln_mesh_pts(:,n);
    for i = 1 : nu
        Y_i_ast(:,i) = y_0;
        for j = 1 : i-2
            fj = odef(t_0 + dT*c_explicit(j),Y_i(:,j));
            Y_i_ast(:,i) = Y_i_ast(:,i)  + dT*A_explicit(i,j)*fj;
        end
        if i > 1
            Y_i_ast(:,i) = Y_i_ast(:,i) + dT*A_explicit(i,i-1)*odef(t_0+dT*c_explicit(i-1)*dT,Y_i(:,i-1));
        end
        
        guess = Y_i_ast(:,i);
        [sol,ithist,ierr] = nsold(guess, @nonlin_sys, [1e-14 1e-14], [40 1 0 1]);
        Y_i(:,i) = sol;
        
        Y_i_all(n,:,i) = Y_i(:,i);
        Y_i_ast_all(n,:,i) = Y_i_ast(:,i);
    end
    
    
    
    
    y_1 = y_0;
    for i = 1 : nu
        y_1 = y_1 +  dT*w_explicit(i)*odef(t_0 + dT*c_explicit(i),Y_i(:,i)) +  dT*w_implicit(i)*odeg(t_0 + dT*c_explicit(i),Y_i(:,i));
    end
    soln_mesh_pts(:,n+1) = y_1;
    
end
    
    
    function ret = nonlin_sys(yi)
        ret = yi - Y_i_ast(:,i);
        
        for jj = 1 : i-1
            gj = odeg(t_0+dT*c_explicit(jj),Y_i(:,jj));
           ret = ret - dT*A_implicit(i,jj)*gj;
        end
        ret = ret - dT*A_implicit(i,i)*odeg(t_0 + dT*c_explicit(i),yi);
    end
    
end