function [Y,t] = strang_splitting_new(funcs,init_vals,h,t_start,t_end)
% Implementation of Strang splitting scheme for solving a system of
% ODEs
% Input:  funcs      - RHSs of the equations as function handles
%         init_vals  - Initial values for the eqns
%         h          - Length of time step
%         t_start    - Initial time
%         t_end      - End time
% Output: Y          - Matrix where each column is the iterates of one eqn
%         t          - Time steps

% Initial values for Y and t
Y = zeros(length(init_vals));

t(1) = t_start;

for i = 1:length(init_vals)
    Y(1,i) = init_vals(i);
end

% Split functions into 2

k = (length(init_vals))/2; % Midpoint index
dx_dt = funcs(1:k);
dp_dt = funcs(k+1:end);

% Number of steps
n = fix((t_end - t_start)/h);

for i = 2:n
    
    % Create intermediate step for x values
    x_intermediate = zeros(1,k);
    for j = 1:length(x_intermediate)
        x_intermediate(j) = Y(i-1,j) + (h/2)*dx_dt{j}(Y(i-1,k+j));
    end
    
    % Fully update p values
    for j = 1:k
        Y(i,k+j) = Y(i-1,k+j) + h*dp_dt{j}(x_intermediate(1),x_intermediate(2),x_intermediate(3),x_intermediate(4),x_intermediate(5),x_intermediate(6),x_intermediate(7),x_intermediate(8),x_intermediate(9),x_intermediate(10),x_intermediate(11),x_intermediate(12),x_intermediate(13),x_intermediate(14),x_intermediate(15),x_intermediate(16),x_intermediate(17),x_intermediate(18),x_intermediate(19),x_intermediate(20),x_intermediate(21),x_intermediate(22),x_intermediate(23),x_intermediate(24),x_intermediate(25),x_intermediate(26),x_intermediate(27),x_intermediate(28),x_intermediate(29),x_intermediate(30));
    end
    
    % Update x values fully
    for j = 1:k
        Y(i,j) = x_intermediate(j) + (h/2)*dx_dt{j}(Y(i,k+j));
    end
    
    t(i) = t(i-1) + h;
    
end