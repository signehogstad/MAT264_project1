function H = hamiltonian_energy(x_p_values,m)
% Function that calculates the energy of a Hamiltonian
% Input: values - Row vector with position and momentum values at a chosen
%                 given time step
%        m      - Masses of the bodies being considered

m_s = 1.98850*10^30;

% Number of bodies being considered (neglects the sun):
% Divides by 6 because only half of the columns are position vectors
% (divide by 2) and each body is associated with 3 position vectors (divide
% by 3
n = length(x_p_values)/6;

% Fetch position values
pos = x_p_values(1:(length(x_p_values)/2));
mom = x_p_values(length(x_p_values)/2+1:end);

% Initial energy value
H = 0;

% Iterating over each body in the system
for j = 1:n
    
    % Starting column in position vector and momentum vector
    k = 3*(j-1) + 1;
    
    % Contribution of momentum for each body
    % First fetch the momentum vector for the jth body
    momentum_vector = mom(j:j+2)';
    
    abs_val_squared = 0;
    for a = 1:length(momentum_vector)
        abs_val_squared = abs_val_squared + (momentum_vector(a))^2;
    end
    
    energy = abs_val_squared/(2*m(j));
    
    H = H + energy;
    
    % Contribution of position
    
    
end


    