function p = momentum(planet_name, position_vectors, masses, dir)
% Function that calculates the dp/dt term for a given planet
% Input:  planet_name      - Name of body being considered (string)
%         position_vectors - Coordinates of all bodies (row vector)
%         masses           - Masses of bodies being considered (vector)
%         dir              - Spatial direction (x=1,y=2,z=3)
% Output: p                - Computed value of dp/dt

% First, determine which body is being considered
if planet_name == 'Earth'
    index = 1;
elseif planet_name == 'Venus'
    index = 2;
elseif planet_name == 'Mrcry'
    index = 3;
elseif planet_name == 'Marss'
    index = 4;
elseif planet_name == 'Jupit'
    index = 5;
elseif planet_name == 'Satur'
    index = 6;
elseif planet_name == 'Uranu'
    index = 7;
elseif planet_name == 'Neptu'
    index = 8;
elseif planet_name == 'Mooon'
    index = 9;
else
    disp('Not a valid body')
end

% Next, fetch positional arguments and mass for the body

pos_main = position_vectors(3*(index-1)+1:3*(index-1)+3);
mass_main = masses(index);

% Remove main positional arguments and mass from input data
position_vectors(3*(index-1)+1:3*(index-1)+3) = [];
masses(index) = [];

% Constants
G = 6.67428*10^-11;
sun_mass = 1.98850*10^30;

% Initialize (P_int = [px py pz])
P_int = zeros(1,3);

% Contribution from the sun in each spatial direction
for j = 1:3
    pos = pos_main;
    sun_contr = -G*sun_mass*mass_main*pos(j)*(pos(1)^2+pos(2)^2+pos(3)^2)^(-3/2);
    P_int(j) = P_int(j) + sun_contr;
end

k = length(position_vectors)/3;
% Contribution from other planets
i = 1;
while i < (k+1)
  
    pos = pos_main;

    % Contribution from other planets
    % Multiply by -G*mass_main as final step
    
    % For each spatial dimensjon (x=1, y=2, z=3)
    for j = 1:3
        contr = masses(i)*(pos(j) - position_vectors(1))*((pos(1)-position_vectors(1))^2+(pos(2)-position_vectors(2))^2+(pos(3)-position_vectors(3))^2)^(-3/2);
        contr = contr * (-G)*mass_main;
        
        P_int(j) = P_int(j) + contr;
        
    end
    
    position_vectors(1:3) =  [];
    i = i + 1;
    
end

P = P_int;

p = P(dir);