function [Y,t,H] = strang_splitting_combined(funcs, init_vals, h, t_start, t_end)
% Implementation of Strang splitting scheme for solving a system of
% ODEs
% Input:  funcs      - RHSs of the equations as function handles
%         init_vals  - Initial values for the eqns
%         h          - Length of time step
%         t_start    - Initial time
%         t_end      - End time
% Output: Y          - Matrix where each column is the iterates of one eqn
%         t          - Time steps

% Masses
m_e = 5.97237*10^24;
m_v = 4.86747*10^24;
m_m = 3.285*10^23;
m_ma = 6.39*10^23;

m = [m_e m_v m_m m_ma];

% Calculating number of time steps
n = ceil((t_end - t_start)/h);

% Initial values for Y and t
Y = zeros(length(init_vals));

for i = 1:length(init_vals)
    Y(1,i) = init_vals(i);
end

t(1) = t_start;

x_p_values = Y(1,:);
H(1) = hamiltonian_energy(x_p_values,m);


for i = 1:n
    
    % Updates time
    t(i+1) = t(i) + h;
    
    % First half of x values (intermediate step)
    
    % Initialize vector of intermediate x values
    x_intermediate = zeros(1,length(init_vals)/2);
    for j = 1:length(x_intermediate)
        k = length(x_intermediate) + j;
        x_intermediate(j) = Y(i,j) + ((t(i+1)-t(i))/2)*funcs{j}(Y(i,k));
    end
    
    % Updating p values
    for j = 1:length(x_intermediate)
        k = length(x_intermediate) + j;
        %Y(i+1,k) = Y(i,k) + (t(i+1)-t(i))*funcs{k}(x_intermediate(1),x_intermediate(2),x_intermediate(3),x_intermediate(4),x_intermediate(5),x_intermediate(6));
        Y(i+1,k) = Y(i,k) + (t(i+1)-t(i))*funcs{k}(x_intermediate(1),x_intermediate(2),x_intermediate(3),x_intermediate(4),x_intermediate(5),x_intermediate(6),x_intermediate(7),x_intermediate(8),x_intermediate(9),x_intermediate(10),x_intermediate(11),x_intermediate(12));
    end


    % Second half of x values
    for j = 1:length(x_intermediate)
        k = length(x_intermediate) + j;
        Y(i+1,j) = x_intermediate(j) + ((t(i+1)-t(i))/2)*funcs{j}(Y(i+1,k));
    end
    
    % Hamiltonian energy
    x_p_values = Y(i+1,:);
    H(i+1) = hamiltonian_energy(x_p_values,m);
    
end