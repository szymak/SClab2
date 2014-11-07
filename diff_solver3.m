function Y = diff_solver3(p0, dt, T, p_fun, method)

% Function diff_solver3 calculate output of the specific 
% differential equation:
%                 p'(t) = fun(p(t))
% using one of three methods:	Explicit Euler,
%                               Method of Heun,
%                               Runge-Kutta method.
% Input:
%                 p0 - initial value p(0),
%                 dt - time step of the approximation,
%                 T - time of the approximation,
%                 p_fun - function p'(t),
%                 method - chosen method:  'euler' - Explicit Euler,
%                                          'heun'  - Method of Heun,
%                                          'r-k'   - Runge-Kutta method.
% Output:
%                 Y - array of calculated values.

    t = 0:dt:T;

    switch (method),    

        case char('euler'),
            % Explicit Euler method
            p_euler = zeros(1,length(t));
            p_euler(1) = p0;
            for i = 1:length(t)-1
                p_euler(i+1) = p_euler(i) + p_fun(p_euler(i))*dt;
            end
            % ----------------------
            Y = p_euler;

        case char('heun'),
            % Method of Heun
            p_heun = zeros(1,length(t));
            p_heun(1) = p0;
            for i = 1:length(t)-1
                p_heun_next = p_heun(i) + p_fun(p_heun(i))*dt;
                p_heun(i+1) = p_heun(i) + (dt/2)*(p_fun(p_heun(i)) + p_fun(p_heun_next));
            end
            % ----------------------
            Y = p_heun;
            
        case char('r-k'),
            %  Runge-Kutta method (fourth order)
            p_rk = zeros(1,length(t));
            p_rk(1) = p0;
            for i=1:length(t)-1
                y1 = p_fun(p_rk(i));
                y2 = p_fun(p_rk(i)+(dt/2)*y1);
                y3 = p_fun(p_rk(i)+(dt/2)*y2);
                y4 = p_fun(p_rk(i)+dt*y3);
                p_rk(i+1) = p_rk(i) + (dt/6)*(y1+2*y2+2*y3+y4);  
            end
            % ----------------------
            Y = p_rk;

    end
    
end