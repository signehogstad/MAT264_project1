function [Y,t] = strang_splitting_3D(funcs, init_vals, h, t_start, t_end)
% Implementation of Strang splitting scheme for solving a system of
% ODEs
% Input:  funcs      - RHSs of the equations as function handles
%         init_vals  - Initial values for the eqns
%         h          - Length of time step
%         t_start    - Initial time
%         t_end      - End time
% Output: Y          - Matrix where each column is the iterates of one eqn
%         t          - Time steps

% Calculating number of time steps
n = ceil((t_end - t_start)/h);

% Initial values for Y and t
Y = zeros(length(init_vals));

for i = 1:length(init_vals)
    Y(1,i) = init_vals(i);
end

t(1) = t_start;


for i = 1:n
    
    % Updates time
    t(i+1) = t(i) + h;
    
    % First half of x values (intermediate step)
    
    % Initialize vector of intermediate x values
    x_intermediate = zeros(1,length(init_vals)/2);
    for j = 1:length(x_intermediate)
        x_intermediate(j) = Y(i,j) + ((t(i+1)-t(i))/2)*funcs{j}(Y(i,1),Y(i,2),Y(i,3),Y(i,4),Y(i,5),Y(i,6),t(i));
    end
    
    % Updating p values
    for j = 1:length(x_intermediate)
        k = length(x_intermediate) + j;
        Y(i+1,k) = Y(i,k) + (t(i+1)-t(i))*funcs{k}(x_intermediate(1),x_intermediate(2),x_intermediate(3),Y(i,4),Y(i,5),Y(i,6),t(i));
    end


    % Second half of x values
    for j = 1:length(x_intermediate)
        Y(i+1,j) = x_intermediate(j) + ((t(i+1)-t(i))/2)*funcs{j}(Y(i,1),Y(i,2),Y(i,3),Y(i+1,4),Y(i+1,5),Y(i+1,6),t(i));
    end
    
end
    
