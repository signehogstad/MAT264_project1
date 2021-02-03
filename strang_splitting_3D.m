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
    x1_intermediate = Y(i,1) + ((t(i+1)-t(i))/2)*funcs{1}(Y(i,1),Y(i,2),Y(i,3),Y(i,4),Y(i,5),Y(i,6),t(i));
    x2_intermediate = Y(i,2) + ((t(i+1)-t(i))/2)*funcs{2}(Y(i,1),Y(i,2),Y(i,3),Y(i,4),Y(i,5),Y(i,6),t(i));
    x3_intermediate = Y(i,3) + ((t(i+1)-t(i))/2)*funcs{3}(Y(i,1),Y(i,2),Y(i,3),Y(i,4),Y(i,6),Y(i,6),t(i));
    
    % Updating p values
    Y(i+1,4) = Y(i,4) + (t(i+1)-t(i))*funcs{4}(x1_intermediate,x2_intermediate,x3_intermediate,Y(i,4),Y(i,5),Y(i,6),t(i));
    Y(i+1,5) = Y(i,5) + (t(i+1)-t(i))*funcs{5}(x1_intermediate,x2_intermediate,x3_intermediate,Y(i,4),Y(i,5),Y(i,6),t(i));
    Y(i+1,6) = Y(i,6) + (t(i+1)-t(i))*funcs{6}(x1_intermediate,x2_intermediate,x3_intermediate,Y(i,4),Y(i,5),Y(i,6),t(i));
    
    % Second half of x values
    Y(i+1,1) = x1_intermediate + ((t(i+1)-t(i))/2)*funcs{1}(Y(i,1),Y(i,2),Y(i,3),Y(i+1,4),Y(i+1,5),Y(i+1,6),t(i));
    Y(i+1,2) = x2_intermediate + ((t(i+1)-t(i))/2)*funcs{2}(Y(i,1),Y(i,2),Y(i,3),Y(i+1,4),Y(i+1,5),Y(i+1,6),t(i));
    Y(i+1,3) = x3_intermediate + ((t(i+1)-t(i))/2)*funcs{3}(Y(i,1),Y(i,2),Y(i,3),Y(i+1,4),Y(i+1,5),Y(i+1,6),t(i));
    
end
    