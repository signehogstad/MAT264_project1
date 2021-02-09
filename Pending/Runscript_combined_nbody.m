hold off
clear
clc

% List of bodies being considered (except the sun)
% 1) Earth
% 2) Venus
% 3) Mercury
% 4) Mars

% Constants
m_s = 1.98847*10^30;
G = 6.67408*10^-11;

% Masses
m_e = 5.97237*10^24;
m_v = 4.8685*10^24;
m_m = 3.3011*10^23;
m_ma = 0.64171*10^24;

m = [m_e, m_v, m_m, m_ma];


% RHSs for displacement and momentum

% Position vectors
x_e = @(px_e) position(px_e,m_e);
y_e = @(py_e) position(py_e,m_e);
z_e = @(pz_e) position(pz_e,m_e);

x_v = @(px_v) position(px_v,m_v);
y_v = @(py_v) position(py_v,m_v);
z_v = @(pz_v) position(pz_v,m_v);

x_m = @(px_m) position(px_m,m_m);
y_m = @(py_m) position(py_m,m_m);
z_m = @(pz_m) position(pz_m,m_m);

x_ma = @(px_ma) position(px_ma,m_ma);
y_ma = @(py_ma) position(py_ma,m_ma);
z_ma = @(pz_ma) position(pz_ma,m_ma);


% Momentum vectors

px_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,1);
py_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,2);
pz_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,3);

px_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,1);
py_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,2);
pz_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,3);

px_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,1);
py_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,2);
pz_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,3);

px_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,1);
py_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,2);
pz_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma],m,3);

funcs = {x_e, y_e, z_e, x_v, y_v, z_v, x_m, y_m, z_m, x_ma, y_ma, z_ma, px_e, py_e, pz_e, px_v, py_v, pz_v, px_m, py_m, pz_m, px_ma, py_ma, pz_ma};


% Initial values

x_init_e = [-2.521092863852298*10^10, 1.449279195712076*10^11, -6.164888475164771*10^5];
x_init_v = [-1.075055502695123*10^11, -3.366520720591562*10^9, 6.159219802771119*10^9];
p_init_e = [-2.983983333368269*m_e*10^4, -5.207633918704476*m_e*10^3, 6.169062303484907*m_e*10^-2];
p_init_v = [8.891598046362434*m_v*10^2, -3.515920774124290*m_v*10^4, -5.318594054684045*m_v*10^2];
x_init_m = [-2.105262111032039*10^10, -6.640663808353403*10^10, -3.492446023382954*10^9];
p_init_m = [3.665298706393840*m_m*10^4, -1.228983810111077*m_m*10^4, -4.368172898981951*m_m*10^3];
x_init_ma = [2.079950549908331*10^11, -3.143009561106971*10^9, -5.178781160069674*10^9];
p_init_ma = [1.295003532851602*m_ma*10^3, 2.629442067068712*m_ma*10^4, 5.190097267545717*m_ma*10^2];


init_vals = [x_init_e, x_init_v, x_init_m, x_init_ma, p_init_e, p_init_v, p_init_m, p_init_ma]; 

h = 50000;

t_start = 0;

t_end = 4000000000;
%t_end = 3721852800; % 2117 Venus transit
%t_end = 392169600; % 2012 Venus transit
%t_end = 139968000; % 2004 Venus transit
%t_end = 31558118; %Earth orbital time
%t_end = 19414166; %Venus orbital time
%t_end = 7600530; % Mercury orbital time

% Try to solve

[Y,t] = strang_splitting_combined_nbody(funcs, init_vals, h, t_start, t_end);


% Make vectors from the sun to the earth and the sun to venus

dir_e = [Y(:,1), Y(:,2), Y(:,3)];

dir_v = [Y(:,4), Y(:,5), Y(:,6)];

% Calculate the inner product

for i = 1:length(Y(:,1))
    dotprod = dot(dir_e(i,:),dir_v(i,:));
    norm1 = norm(dir_e(i,:));
    norm2 = norm(dir_v(i,:));
    theta(i) = acos(dotprod/(norm1*norm2));
end

theta(end);
k = 1;

for i = 1:length(theta)
    if theta(i) < 5*10^-3
        index(k) = i;
        theta_small(k) = theta(i);
        k = k + 1;
    end
end

% error_x = abs(Y(1,7) - Y(end,7));
% error_y = abs(Y(1,8) - Y(end,8));
% error_z = abs(Y(1,9) - Y(end,9));
% feilvektor = [error_x error_y error_z]';
% 
% abs_feil = sqrt(error_x^2 + error_y^2 + error_z^2);

% Plotting

%  h = animatedline;
%  d = animatedline;
%  e = animatedline;
%  w = animatedline;
%  m_1 = animatedline; % mercury
%  m_2 = animatedline; % mercury
%  ma_1 = animatedline; % mars
%  ma_2 = animatedline; % mars
%  axis([-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12])
%  axis square
%  
%  hold on
%  plot(0,0,'r*','Markersize',5)
%  text(0,0,'sun')
%  
%  
%  for k = 1:length(Y(:,1))
%      addpoints(h,Y(k,1),Y(k,2),Y(k,3));
%      addpoints(d,Y(k,4),Y(k,5),Y(k,6));
%      addpoints(m_1,Y(k,7),Y(k,8),Y(k,9));
%      addpoints(ma_1,Y(k,10),Y(k,11),Y(k,12));
%      set(e,'Color','r','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
%      addpoints(e,Y(k,1),Y(k,2),Y(k,3));
%      set(w,'Color','b','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
%      addpoints(w,Y(k,4),Y(k,5),Y(k,6));
%      set(m_2,'Color','k','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
%      addpoints(m_2,Y(k,7),Y(k,8),Y(k,9));
%      set(ma_2,'Color','g','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
%      addpoints(ma_2,Y(k,10),Y(k,11),Y(k,12));
%      drawnow
%  end
% 
% H = 0;
% for i = 1:length(Y(:,1))
%     H(i) = hamiltonian_energy(Y(i,:),m);
% end
