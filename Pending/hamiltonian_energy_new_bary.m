function H = hamiltonian_energy_new_bary(x_p_val,m)
% Function that calculates the energy of a Hamiltonian system
% Input:  x_p_val - Row vector with position and momentum values at a chosen
%                   given time step
%         m       - Masses of the bodies being considered
% Output: H       - Hamiltonian at the given time step

G = 6.67428*10^-11;

% Number of bodies being considered
n = length(x_p_val)/6; % Each body is associated w/ 6 vectors

% Fetch position values and momentum values
pos = x_p_val(1:(length(x_p_val)/2));
mom = x_p_val(length(x_p_val)/2+1:end);

% Initial energy value
H  = 0;

% Iterating over each body in the system
for j = 1:n
    
    % First calculate momentum contribution for each body
    
    % Starting column in position vector and momentum vector
    k = 3*(j-1) + 1;
    
    % Momentum vector for body i
    mom_vec = mom(k:k+2);
    
    % Expression for momentum contribution
    numerator = mom_vec(1)^2 + mom_vec(2)^2 + mom_vec(3)^2;
    denominator = 2*m(j);
    
    mom_contr = numerator/denominator;
    
    H = H + mom_contr;
    
    % Next, calculate gravity contribution for each body
    
    m_copy = m;
    mass_main = m(j); % Fetch mass of current body
    m_copy(j) = 0; % Set to 0 in mass vector so l=j term will cancel
    
    grav_contr = 0;
    for i = 1:n
        mass_product = mass_main*m_copy(i);
        if mass_product == 0
            contr = 0;
        else
            dist_vec = [pos(k) - pos(3*(i-1) + 1), pos(k+1) - pos(3*(i-1) + 2), pos(k+2) - pos(3*(i-1) + 3)];
            abs_val = norm(dist_vec);
            contr = mass_product/abs_val;
        end
        grav_contr = grav_contr + contr;
    end
    
    grav_contr = grav_contr*(-G);
 
    H = H + grav_contr;
    
end
    