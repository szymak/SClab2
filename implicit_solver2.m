function Y = implicit_solver2(p0, dt, T, eps_min, nmax, fun_newton, fun_newton_prim, method)

% Function implicit_solver2 calculate output of the specific 
% differential equation:
%                 p'(t) = fun(p(t))
% using one of two methods:	Implicit Euler,
%                           Adams-Moulton method.
% Input:
%                 p0 - initial value,
%                 dt - timestep,
%                 T - T end,
%                 eps_min - accuracy in Newton's Method,
%                 nmax - maximal number of steps for Newton's method,
%                 fun_newton - function needed for Newton's method,
%                 fun_newton_prim - derivative of the above function "p_fun_newton" also 
%                                   needed for Newton's method,
%                 method - name of used method.
% Output:
%                 Y - array of calculated values.

    t = 0:dt:T;
    newton_0 = 20; % start point of Newton's method
    
    switch (method),    

        case char('euler')
            % Implicit Euler method
            p_euler = zeros(1,length(t));
            p_euler(1) = p0;
            for i=1:length(t)-1,
                [p_solution, if_found] = newton_solver(newton_0,p_euler(i),dt,eps_min,nmax,fun_newton,fun_newton_prim);            
                if if_found == 0 % if_found determines whether solution has been found
                    p_euler = NaN(1,length(t));
                    break;
                end
                p_euler(i+1) = p_solution;
                newton_0 = p_solution;
            end
            % ----------------------
            Y = p_euler;

        case char('adams')
            % Adams-Moulton method
            p_adams = zeros(1,length(t));
            p_adams(1) = p0;
            for i=1:length(t)-1,
                [p_solution, if_found] = newton_solver(newton_0,p_adams(i),dt,eps_min,nmax,fun_newton,fun_newton_prim);
                if if_found == 0 % if_found determines whether solution has been found
                    p_adams = NaN(1,length(t));
                    break;
                end  
                p_adams(i+1) = p_solution;
                newton_0 = p_solution;
            end
            % ----------------------
            Y = p_adams;

    end

end