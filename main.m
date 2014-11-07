%       Scientific Computing lab
%             Worksheet 2      

clc; clear all; close all

% changeable part --------------------------------------

p0 = 20; % initial function value
T = 5; % time of simulation
dt_tab = [1/2 1/4 1/8 1/16 1/32]; % different time steps
eps_min = 1e-4; % accuracy in Newton's Method
nmax = 100; % maximal number of steps for Newton's method

p_fun = @(p)(7*(1-(p/10))*p); % basic equation
p_analytical = @(tp)(200./(20-10*exp(-7*tp))); % analytical solution

% function needed for Newton's method (Euler's method)
p_fun_newton = @(p_nplus1,p_n,dt)(0.7*dt*p_nplus1.*p_nplus1+(1-7*dt)*p_nplus1-p_n); 
% derivative of the above function "p_fun_newton" also needed for Newton's method
p_fun_newton_prim = @(p_nplus1,dt)(1.4*dt*p_nplus1+1-7*dt); 

% function needed for Newton's method (Adams-Moulton Method method)
p_fun_newton2 = @(p_nplus1,p_n,dt)(0.35*dt*p_nplus1.*p_nplus1+(1-3.5*dt)*p_nplus1+0.35*dt*p_n*p_n - (1+3.5*dt)*p_n); 
% derivative of the above function "p_fun_newton2" also needed for Newton's method
p_fun_newton_prim2 = @(p_nplus1,dt)(0.7*dt*p_nplus1+1-3.5*dt); 

% 1st linearisation function expressed as y(n+1)=f(y(n))
p_fun_lin1 = @(y,dt)((-0.35*dt*y*y+(1+7*dt)*y)/(1+0.35*dt*y));
% 2nd linearisation function expressed as y(n+1)=f(y(n))
p_fun_lin2 = @(y,dt)((-0.35*dt*y*y+(1+3.5*dt)*y)/(1-3.5*dt+0.35*dt*y));

% -----------------------------------------------------

E_fun = @(dt,T,p_num,p_exact)(sqrt((dt/T).*sum((p_num-p_exact).^2))); % function for calculating approximation error
t_graph = 0:0.01:T; % time with very small time step just for the nice graph of analytical solution
p_analytical_graph = p_analytical(t_graph); % analytical solution graph data
row_names = {'tstep', 'error', 'error red.'};
variable_names = cellfun( @(i) ['ts' num2str(i)], num2cell(1:length(dt_tab)), 'UniformOutput', false);

dt_tab = sort(dt_tab,'descend'); % sorting elements to ensure that the smallest step is in the end

% Explicit methods ----------

methods = struct('name', {'Explicit Euler', 'Method of Heun'}, 'id', {'euler', 'heun'});

for method = methods
    
    data = zeros(3, length(dt_tab)); % matrix with results
    
    figure('name', method.name);
    
    for i = length(dt_tab):-1:1

        dt = dt_tab(i);
        t = 0:dt:T; 

        Y = diff_solver3(p0, dt, T, p_fun, method.id);

        E1 = E_fun(dt, T, Y, p_analytical(t)); % first error

        data(:,i) = [dt; E1; 1];

        % plotting
        subplot(ceil(length(dt_tab)/2), 2, i);
        plot(t_graph, p_analytical_graph, 'b-', t, Y, 'r');
        title(['dt = ' num2str(dt_tab(i))],'fontweight','bold');
        legend('Analytical solution', method.name);        
        grid on; ylabel('y'); xlabel('t');  
        axis([0 T -5 20]);
        
    end

    % error reducing factor
    for i = 2:length(dt_tab)
        data(3,i) = data(2,i-1) / data(2,i);
    end

    % visualize table
    data = array2table(data, 'RowNames', row_names, 'VariableNames', variable_names);
    disp(method.name); disp(data);
    
end

% --------------------------------

% Implicit methods ---------------

methods2 = struct('name', {'Implicit Euler', 'Adams-Moulton method'}, 'id', {'euler', 'adams'});

for method=methods2,

    data = zeros(3, length(dt_tab)); % matrix with results
    
    figure('name', method.name);    
    
    for i = 1:length(dt_tab)

        dt = dt_tab(i);
        t = 0:dt:T; 
        if strcmp(method.id,'euler'),
            Y = implicit_solver2(p0, dt, T, eps_min,nmax, p_fun_newton, p_fun_newton_prim, method.id);
        else
            Y = implicit_solver2(p0, dt, T, eps_min,nmax, p_fun_newton2, p_fun_newton_prim2, method.id);
        end
        
        E1 = E_fun(dt, T, Y, p_analytical(t)); % first error
        data(:,i) = [dt; E1; 1];
        
        % plotting
        subplot(ceil(length(dt_tab)/2), 2, i);
        plot(t_graph, p_analytical_graph, 'b-', t, Y, 'r');
        title(['dt = ' num2str(dt_tab(i))],'fontweight','bold');
        legend('Analytical solution', method.name);        
        grid on; ylabel('y'); xlabel('t');  
        axis([0 T -5 20]);
        
    end
    
    % error reducing factor
    for i = 2:length(dt_tab)
        data(3,i) = data(2,i-1) / data(2,i);
    end
    
    % visualize table
    data = array2table(data, 'RowNames', row_names, 'VariableNames', variable_names);
    disp(method.name); disp(data);

end

% ---------------------------------

% Linearisation methods -----------

methods3 = struct('name', {'Adams-Moulton 1st linearisation', 'Adams-Moulton 2nd linearisation'},...
                                                                    'id', {'adams1', 'adams2'});
                                                                
for method=methods3,
    
    data = zeros(3, length(dt_tab)); % matrix with results
    
    figure('name', method.name);
        
    for i = 1:length(dt_tab)

        dt = dt_tab(i);
        t = 0:dt:T; 
        Y_lin = zeros(1,length(t));
        Y_lin(1) = p0;

        for k=1:length(t)-1,
            if strcmp(method.id,'adams1'),
                Y_lin(k+1) = p_fun_lin1(Y_lin(k),dt);
            else
                Y_lin(k+1) = p_fun_lin2(Y_lin(k),dt);
            end
        end
        E1 = E_fun(dt, T, Y_lin, p_analytical(t)); % first error
        data(:,i) = [dt; E1; 1];

        subplot(ceil(length(dt_tab)/2), 2, i);
        plot(t_graph, p_analytical_graph, 'b-', t, Y_lin, 'r');
        title(['dt = ' num2str(dt_tab(i))],'fontweight','bold');
        legend('Analytical solution', method.name);        
        grid on; ylabel('y'); xlabel('t');  
        axis([0 T -5 20]);

    end
        
    % error reducing factor
    for i = 2:length(dt_tab)
        data(3,i) = data(2,i-1) / data(2,i);
    end
    
    % visualize table
    data = array2table(data, 'RowNames', row_names, 'VariableNames', variable_names);
    disp(method.name); disp(data);
    
end

% ---------------------------------------

% Stable cases
row_names2 = {'dt = 1/2', 'dt = 1/4', 'dt = 1/8', 'dt = 1/16', 'dt = 1/32'};
variable_names2 = {'explicit_Euler', 'Heun', 'implicit_Euler', ...
                   'Adams_Moulton', 'Adams_Moulton_l1', 'Adams_Moulton_l2'};
               
stability_data = [' ' ' ' 'x' ' ' 'x' ' ';...
                  ' ' ' ' 'x' 'x' 'x' 'x';...
                  'x' 'x' 'x' 'x' 'x' 'x';...
                  'x' 'x' 'x' 'x' 'x' 'x';...
                  'x' 'x' 'x' 'x' 'x' 'x'];

stability_data = array2table(stability_data, 'RowNames', row_names2, 'VariableNames', variable_names2);
disp('Stable cases'); 
disp(stability_data);
disp('If the tables do not look pretty, change font in command window to "Monospaced".')