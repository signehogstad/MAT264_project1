function [p] = momentum(position,M,m_r,G)
% Function that calculates the momentum in x, y and z direction based on
% the position of other planets and the planet itself, the masses and the
% gravitational constant
p = 0;
for m in m_r
    p = p -G*M*m*position{
p = -G*m(1)*(m_s*x_e*(x_e^2+y_e^2+z_e^2)^(-3/2) + m(2)*(x_e-x_v)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));
end