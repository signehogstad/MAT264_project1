function [Y,t,H] = euler_forward(funcs, init_vals, h, t_start, t_end)
% Implementation of forward (explicit) Euler method for solving a system of
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

%position_funcs = funcs{1:length(init_vals)/2};
%momentum_funcs = funcs{length(init_vals)/2:end};

for i = 1:n
    % Updating position functions
    for k = 1:length(init_vals)/2
        Y(i+1,k) = Y(i,k) + h*funcs{k}(Y(i,length(init_vals)/2+k));
    end
    
    % Updating momentum functions
    for k = (length(init_vals)/2+1):length(init_vals)
        Y(i+1,k) = Y(i,k) + h*funcs{k}(Y(i,1),Y(i,2),Y(i,3),Y(i,4),Y(i,5),Y(i,6),Y(i,7),Y(i,8),Y(i,9),Y(i,10),Y(i,11),Y(i,12));
    end
    
    
    % Updating y for each function
%     Y(i+1,1) = Y(i,1) + h*funcs{1}(Y(i,1),Y(i,2),Y(i,3),Y(i,4));
%     Y(i+1,2) = Y(i,2) + h*funcs{2}(Y(i,1),Y(i,2),Y(i,3),Y(i,4));
%     Y(i+1,3) = Y(i,3) + h*funcs{3}(Y(i,1),Y(i,2),Y(i,3),Y(i,4));
%     Y(i+1,4) = Y(i,4) + h*funcs{4}(Y(i,1),Y(i,2),Y(i,3),Y(i,4));
    % Add more lines + arguments if there are more than 2 functions
    
    % Updating t
    t(i+1) =  t(i) + h;
    
    % Hamiltonian energy
    x_p_values = Y(i+1,:);
    H(i+1) = hamiltonian_energy(x_p_values,m);
    
end