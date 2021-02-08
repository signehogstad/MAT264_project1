hold off
clear
clc

% List of bodies being considered (except the sun)
% 1) Earth
% 2) Venus
% 3) Mercury
% 4) Mars


% Constants
m_s = 1.98850*10^30;
G = 6.67428*10^-11;

% Masses
m_e = 5.97237*10^24;
m_v = 4.86747*10^24;
m_m = 3.285*10^23;
m_ma = 6.39*10^23;

m = [m_e m_v m_m m_ma];


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

x_init_e = [-2.627892928682480*10^10, 1.445102393586391*10^11, 3.022818135935813*10^7];
x_init_v = [-1.085735509178141*10^11, -3.784200933160055*10^9, 6.190064472977990*10^9];
p_init_e = [-2.983052803283506*m_e*10^4, -5.220465685407924*m_e*10^3, -1.014621798034465*m_e*10^-1];
p_init_v = [8.984651054838754*m_v*10^2, -3.517203950794635*m_v*10^4, -5.320225582712421*m_v*10^2];
x_init_m = [-2.212062175862221*10^10, -6.682431829610253*10^10, -3.461601353176080*10^9];
p_init_m = [3.666229236478603*m_m*10^4, -1.230266986781422*m_m*10^4, -4.368336051784789*m_m*10^3];
x_init_ma = [2.069270543147017*10^11, -3.560689745239088*10^9, -5.147936537447235*10^9];
p_init_ma = [1.304308833322233*m_ma*10^3, 2.628158890420931*m_ma*10^4, 5.188465740839767*m_ma*10^2];


init_vals = [x_init_e, x_init_v, x_init_m, x_init_ma, p_init_e, p_init_v, p_init_m, p_init_ma]; 

h = 50000;

t_start = 0;

t_end = 450000000;

% Try to solve

[Y,t,H] = strang_splitting_combined(funcs, init_vals, h, t_start, t_end);


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


% Plotting

 h = animatedline;
 d = animatedline;
 e = animatedline;
 w = animatedline;
 m_1 = animatedline; % mercury
 m_2 = animatedline; % mercury
 ma_1 = animatedline; % mars
 ma_2 = animatedline; % mars
 axis([-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12])
 axis square
 
 hold on
 plot(0,0,'r*','Markersize',5)
 text(0,0,'sun')
 
 
 for k = 1:length(Y(:,1))
     addpoints(h,Y(k,1),Y(k,2),Y(k,3));
     addpoints(d,Y(k,4),Y(k,5),Y(k,6));
     addpoints(m_1,Y(k,7),Y(k,8),Y(k,9));
     addpoints(ma_1,Y(k,10),Y(k,11),Y(k,12));
     set(e,'Color','r','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(e,Y(k,1),Y(k,2),Y(k,3));
     set(w,'Color','b','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(w,Y(k,4),Y(k,5),Y(k,6));
     set(m_2,'Color','k','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(m_2,Y(k,7),Y(k,8),Y(k,9));
     set(ma_2,'Color','g','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(ma_2,Y(k,10),Y(k,11),Y(k,12));
     drawnow
 end
