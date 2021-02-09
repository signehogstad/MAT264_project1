hold off
clear
clc

% List of bodies being considered (except the sun)
% 1) Earth
% 2) Venus
% 3) Mercury
% 4) Mars
% 5) Jupiter
% 6) Saturn
% 7) Uranus
% 8) Neptun

% Constants
m_S = 1.98847*10^30;
G = 6.67408*10^-11;

% Masses
m_e = 5.97237*10^24;
m_v = 4.86747*10^24;
m_m = 3.285*10^23;
m_ma = 6.39*10^23;
m_j = 1898.187*10^24;
m_s = 568.3174*10^24;
m_u = 86.8127*10^24;
m_n = 102.4126*10^24;

m = [m_e m_v m_m m_ma m_j m_s m_u m_n];


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

x_j = @(px_j) position(px_j,m_j);
y_j = @(py_j) position(py_j,m_j);
z_j = @(pz_j) position(pz_j,m_j);

x_s = @(px_s) position(px_s,m_s);
y_s = @(py_s) position(py_s,m_s);
z_s = @(pz_s) position(pz_s,m_s);

x_u = @(px_u) position(px_u,m_u);
y_u = @(py_u) position(py_u,m_u);
z_u = @(pz_u) position(pz_u,m_u);

x_n = @(px_n) position(px_n,m_n);
y_n = @(py_n) position(py_n,m_n);
z_n = @(pz_n) position(pz_n,m_n);


% Momentum vectors

px_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

px_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,1);
py_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,2);
pz_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n) momentum('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n],m,3);

funcs = {x_e, y_e, z_e, x_v, y_v, z_v, x_m, y_m, z_m, x_ma, y_ma, z_ma, x_j, y_j, z_j, x_s, y_s, z_s, x_u, y_u, z_u, x_n, y_n, z_n, px_e, py_e, pz_e, px_v, py_v, pz_v, px_m, py_m, pz_m, px_ma, py_ma, pz_ma, px_j, py_j, pz_j, px_s, py_s, pz_s, px_u, py_u, pz_u, px_n, py_n, pz_n};

% Initial values

x_init_e = [-2.627892928682480*10^10, 1.445102393586391*10^11, 3.022818135935813*10^7];
p_init_e = [-2.983052803283506*m_e*10^4, -5.220465685407924*m_e*10^3, -1.014621798034465*m_e*10^-1];
x_init_v = [-1.085735509178141*10^11, -3.784200933160055*10^9, 6.190064472977990*10^9];
p_init_v = [8.984651054838754*m_v*10^2, -3.517203950794635*m_v*10^4, -5.320225582712421*m_v*10^2];
x_init_m = [-2.212062175862221*10^10, -6.682431829610253*10^10, -3.461601353176080*10^9];
p_init_m = [3.666229236478603*m_m*10^4, -1.230266986781422*m_m*10^4, -4.368336051784789*m_m*10^3];
x_init_ma = [2.069270543147017*10^11, -3.560689745239088*10^9, -5.147936537447235*10^9];
p_init_ma = [1.304308833322233*m_ma*10^3, 2.628158890420931*m_ma*10^4, 5.188465740839767*m_ma*10^2];
x_init_j = [5.978411018921824E+11, 4.387049769612203E+11, -1.520170016953275E+10];
p_init_j = [-m_j*7.892151874487445E+03, m_j*1.114945115750926E+04, m_j*1.304862310943014E+02];
x_init_s = [9.576381404895931E+11, 9.821477386679871E+11, -5.518990106544006E+10];
p_init_s = [-m_s*7.420583259777167E+03, m_s*6.725272210562052E+03, m_s*1.779635671891335E+02];
x_init_u = [2.157706693154863E+12, -2.055242913196469E+12, -3.559266720965648E+10];
p_init_u = [m_u*4.646804486524291E+03, m_u*4.614387135494279E+03, -m_u*4.306254213333638E+01];
x_init_n = [2.513785356833091E+12, -3.739265096828156E+12, 1.907035719387865E+10];
p_init_n = [m_n*4.475371931235400E+03, m_n*3.063567701145199E+03, -m_n*1.662252948584049E+02];

init_vals = [x_init_e, x_init_v, x_init_m, x_init_ma, x_init_j, x_init_s, x_init_u, x_init_n, p_init_e, p_init_v, p_init_m, p_init_ma, p_init_j, p_init_s, p_init_u, p_init_n]; 

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

h1 = animatedline; % earth line
 d = animatedline; % venus line
 m_1 = animatedline; % mercury line
 ma_1 = animatedline; % mars line
 q = animatedline; % jupiter line
 r = animatedline; % saturn line
 s = animatedline; % uranus line
 t1 = animatedline; % neptun line
 
 e = animatedline; % earth dot
 w = animatedline; % venus dot
 m_2 = animatedline; % mercury dot
 ma_2 = animatedline; % mars dot
 q2 = animatedline; % jupiter dot
 r2 = animatedline; % saturn dot
 s2 = animatedline; % uranus dot
 t2 = animatedline; % neptun dot
 
%  axis([-0.25*10^13,0.25*10^13,-0.25*10^13,0.25*10^13,-0.25*10^13,0.25*10^13])
%  axis square
%  
%  hold on
%  plot(0,0,'r*','Markersize',5)
%  text(0,0,'sun')
%  
%  
%  for k = 1:length(Y(:,1))
%      addpoints(h1,Y(k,1),Y(k,2),Y(k,3));
     addpoints(d,Y(k,4),Y(k,5),Y(k,6));
     addpoints(m_1,Y(k,7),Y(k,8),Y(k,9));
     addpoints(ma_1,Y(k,10),Y(k,11),Y(k,12));
     addpoints(q,Y(k,13),Y(k,14),Y(k,15));
     addpoints(r,Y(k,16),Y(k,17),Y(k,18));
     addpoints(s,Y(k,19),Y(k,20),Y(k,21));
     addpoints(t1,Y(k,22),Y(k,23),Y(k,24));
     set(e,'Color','r','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(e,Y(k,1),Y(k,2),Y(k,3));
     set(w,'Color','b','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(w,Y(k,4),Y(k,5),Y(k,6));
     set(m_2,'Color','k','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(m_2,Y(k,7),Y(k,8),Y(k,9));
     set(ma_2,'Color','g','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(ma_2,Y(k,10),Y(k,11),Y(k,12));
     set(q2,'Color','c','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(q2,Y(k,13),Y(k,14),Y(k,15));
     set(r2,'Color','m','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(r2,Y(k,16),Y(k,17),Y(k,18));
     set(s2,'Color','y','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(s2,Y(k,19),Y(k,20),Y(k,21));
     set(t2,'Color','r','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(t2,Y(k,22),Y(k,23),Y(k,24));
     drawnow
%  end
% 
% H = 0;
% for i = 1:length(Y(:,1))
%     H(i) = hamiltonian_energy(Y(i,:),m);
% end
