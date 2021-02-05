hold off
clear
clc


% Constants
m_e = 5.97237*10^24;
m_v = 4.86747*10^24;
m_s = 1.98850*10^30;
G = 6.67428*10^-11;
r = 1.621286506*10^11;


% RHSs for displacement and momentum

% Position vectors
x_e = @(px_e) px_e*(1/m_e);
y_e = @(py_e) py_e*(1/m_e);
z_e = @(pz_e) pz_e*(1/m_e);

x_v = @(px_v) px_v*(1/m_v);
y_v = @(py_v) py_v*(1/m_v);
z_v = @(pz_v) pz_v*(1/m_v);

% Momentum vectors
%px_e = @(x_e,y_e,z_e,x_v,y_v,z_v) ((-G*m_e*m_s)*(x_e^2+y_e^2+z_e^2)^(-3/2)*(x_e));
%py_e = @(x_e,y_e,z_e,x_v,y_v,z_v) ((-G*m_e*m_s)*(x_e^2+y_e^2+z_e^2)^(-3/2)*(y_e));
%pz_e = @(x_e,y_e,z_e,x_v,y_v,z_v) ((-G*m_e*m_s)*(x_e^2+y_e^2+z_e^2)^(-3/2)*(z_e));

px_e = @(x_e,y_e,z_e,x_v,y_v,z_v) -G*m_e*(m_s*x_e*(x_e^2+y_e^2+z_e^2)^(-3/2) + m_v*(x_e-x_v)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));
py_e = @(x_e,y_e,z_e,x_v,y_v,z_v) -G*m_e*(m_s*y_e*(x_e^2+y_e^2+z_e^2)^(-3/2) + m_v*(y_e-y_v)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));
pz_e = @(x_e,y_e,z_e,x_v,y_v,z_v) -G*m_e*(m_s*z_e*(x_e^2+y_e^2+z_e^2)^(-3/2) + m_v*(z_e-z_v)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));

px_v = @(x_e,y_e,z_e,x_v,y_v,z_v) -G*m_v*(m_s*x_v*(x_v^2+y_v^2+z_v^2)^(-3/2) + m_e*(x_v-x_e)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));
py_v = @(x_e,y_e,z_e,x_v,y_v,z_v) -G*m_v*(m_s*y_v*(x_v^2+y_v^2+z_v^2)^(-3/2) + m_e*(y_v-y_e)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));
pz_v = @(x_e,y_e,z_e,x_v,y_v,z_v) -G*m_v*(m_s*z_v*(x_v^2+y_v^2+z_v^2)^(-3/2) + m_e*(z_v-z_e)*((x_e-x_v)^2+(y_e-y_v)^2+(z_e-z_v)^2)^(-3/2));

%px_v = @(x_e,y_e,z_e,x_v,y_v,z_v) ((-G*m_v*m_s)*(x_v^2+y_v^2+z_v^2)^(-3/2)*(x_v));
%py_v = @(x_e,y_e,z_e,x_v,y_v,z_v) ((-G*m_v*m_s)*(x_v^2+y_v^2+z_v^2)^(-3/2)*(y_v));
%pz_v = @(x_e,y_e,z_e,x_v,y_v,z_v) ((-G*m_v*m_s)*(x_v^2+y_v^2+z_v^2)^(-3/2)*(z_v));


funcs = {x_e, y_e, z_e, x_v, y_v, z_v, px_e, py_e, pz_e, px_v, py_v, pz_v};


% Initial values

init_vals = [-2.627892928682480*10^10, 1.445102393586391*10^11, 3.022818135935813*10^7, -1.085735509178141*10^11, -3.784200933160055*10^9, 6.190064472977990*10^9,-2.983052803283506*m_e*10^4, -5.220465685407924*m_e*10^3, -1.014621798034465*m_e*10^-1, 8.984651054838754*m_v*10^2, -3.517203950794635*m_v*10^4, -5.320225582712421*m_v*10^2]; 

h = 10000;

t_start = 0;

t_end = 400000000;

% Try to solve

[Y,t] = strang_splitting_combined(funcs, init_vals, h, t_start, t_end);


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

 h = animatedline;
 d = animatedline;
 e = animatedline;
 w = animatedline;
 axis([-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12])
 axis square
 
 hold on
 plot(0,0,'r*','Markersize',5)
 text(0,0,'sun')
 
 
 for k = 1:length(Y(:,1))
     addpoints(h,Y(k,1),Y(k,2),Y(k,3));
     addpoints(d,Y(k,4),Y(k,5),Y(k,6));
     set(e,'Color','r','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(e,Y(k,1),Y(k,2),Y(k,3));
     set(w,'Color','b','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
     addpoints(w,Y(k,4),Y(k,5),Y(k,6));
     drawnow
 end
