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

% Earth
fx1 = @(x1,x2,x3,p1,p2,p3,t) p1*(1/m_e);
fx2 = @(x1,x2,x3,p1,p2,p3,t) p2*(1/m_e);
fx3 = @(x1,x2,x3,p1,p2,p3,t) p3*(1/m_e);
fp1 = @(x1,x2,x3,p1,p2,p3,t) ((-G*m_e*m_s)*(x1^2+x2^2+x3^2)^(-3/2)*(x1));
fp2 = @(x1,x2,x3,p1,p2,p3,t) ((-G*m_e*m_s)*(x1^2+x2^2+x3^2)^(-3/2)*(x2));
fp3 = @(x1,x2,x3,p1,p2,p3,t) ((-G*m_e*m_s)*(x1^2+x2^2+x3^2)^(-3/2)*(x3));

% Venus
fx1_v = @(x1,x2,x3,p1,p2,p3,t) p1*(1/m_v);
fx2_v = @(x1,x2,x3,p1,p2,p3,t) p2*(1/m_v);
fx3_v = @(x1,x2,x3,p1,p2,p3,t) p3*(1/m_v);
fp1_v = @(x1,x2,x3,p1,p2,p3,t) ((-G*m_v*m_s)*(x1^2+x2^2+x3^2)^(-3/2)*(x1));
fp2_v = @(x1,x2,x3,p1,p2,p3,t) ((-G*m_v*m_s)*(x1^2+x2^2+x3^2)^(-3/2)*(x2));
fp3_v = @(x1,x2,x3,p1,p2,p3,t) ((-G*m_v*m_s)*(x1^2+x2^2+x3^2)^(-3/2)*(x3));


funcs = {fx1, fx2, fx3, fp1, fp2, fp3};

funcs_v = {fx1_v, fx2_v, fx3_v, fp1_v, fp2_v, fp3_v};


% Initial values

init_vals = [-2.627892928682480*10^10, 1.445102393586391*10^11, 3.022818135935813*10^7, -2.983052803283506*m_e*10^4, -5.220465685407924*m_e*10^3, -1.014621798034465*m_e*10^-1]; 
init_vals_v = [-1.085735509178141*10^11, -3.784200933160055*10^9, 6.190064472977990*10^9, 8.984651054838754*m_v*10^2, -3.517203950794635*m_v*10^4, -5.320225582712421*m_v*10^2];

h = 1000000;

t_start = 0;

t_end = 2000000000;

% Try to solve

[Y,t] = strang_splitting_3D(funcs, init_vals, h, t_start, t_end);

[Y_v,t_v] = strang_splitting_3D(funcs_v, init_vals_v, h, t_start, t_end);

% Make vectors from the sun to the earth and the sun to venus

dir_e = [Y(:,1), Y(:,2), Y(:,3)];

dir_v = [Y_v(:,1), Y_v(:,2), Y_v(:,3)];

% Calculate the inner product

for i = 1:length(Y(:,1))
    dotprod = dot(dir_e(i,:),dir_v(i,:));
    norm1 = norm(dir_e(i,:));
    norm2 = norm(dir_v(i,:));
    theta(i) = acos(dotprod/(norm1*norm2));
end

h = animatedline;
h_v = animatedline;
d = animatedline;
d_v = animatedline;
axis([-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12,-0.25*10^12,0.25*10^12])
axis square

hold on
plot(0,0,'r*','Markersize',5)
text(0,0,'sun')


for k = 1:length(Y(:,1))
    addpoints(h,Y(k,1),Y(k,2),Y(k,3));
    addpoints(h_v,Y_v(k,1),Y_v(k,2),Y_v(k,3));
    set(d,'Color','r','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
    addpoints(d,Y(k,1),Y(k,2),Y(k,3));
    set(d_v,'Color','b','Marker','.','MarkerSize',25,'MaximumNumPoints',1);
    addpoints(d_v,Y_v(k,1),Y_v(k,2),Y_v(k,3));
    drawnow
    %hold on
    %plot(Y(k,1),Y(k,2),'or','MarkerSize',2,'MarkerFaceColor','r')
    %plot(Y_v(k,1),Y_v(k,2),'or','MarkerSize',2,'MarkerFaceColor','b')
end


%axis square
%hold on
%plot(Y(k,1),Y(k,2),'or','MarkerSize',2,'MarkerFaceColor','r')
