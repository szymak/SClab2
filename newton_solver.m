function [x_zero, if_found] = newton_solver(x0, x_moment, dt, eps_min, nmax, fun, fun_prim)

% Newton's Method
% x0 - initial starting point of Newton's method,
% x_moment - current value of the function in specific time step, 
% dt - timestep,
% eps_min - minimal value close to 0 which will be sufficient,
% nmax - maximal number of steps of the method,
% fun - considered function,
% fun_prim - considered 1st derivative of the function.

% x_zero - calculated root of function,
% if_found - value "0" means that within "nmax" steps the estimated
%            solution was not good enough ( >eps_min )

    x = x0; %Initial 0 point
    n = 0;  %Number of steps
    
    while abs(fun(x,x_moment,dt)) > eps_min && n < nmax
        x = x - fun(x,x_moment,dt) / fun_prim(x,dt);
        n = n+1;
    end
    
    if fun(x,x_moment,dt) > eps_min
        x_zero = 0;
        if_found = 0;
    else
        x_zero = x;
        if_found = 1;
    end

end