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
% 9) Moon

% Constants
m_S = 1.98847*10^30;
G = 6.67408*10^-11;

% Masses
m_e = 5.97219*10^24;
m_v = 48.685*10^23;
m_m = 3.302*10^23;
m_ma = 6.39*10^23;
m_j = 1898.187*10^24;
m_s = 568.3174*10^24;
m_u = 86.8127*10^24;
m_n = 102.4126*10^24;
m_mo = 7.34767309*10^22;
%m_mo = 1;

m = [m_e m_v m_m m_ma m_j m_s m_u m_n m_mo];


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

x_mo = @(px_mo) position(px_mo,m_mo);
y_mo = @(py_mo) position(py_mo,m_mo);
z_mo = @(pz_mo) position(pz_mo,m_mo);


% Momentum vectors

px_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

px_mo = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Mooon',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,1);
py_mo = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Mooon',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,2);
pz_mo = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo) momentum('Mooon',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo],m,3);

funcs = {x_e, y_e, z_e, x_v, y_v, z_v, x_m, y_m, z_m, x_ma, y_ma, z_ma, x_j, y_j, z_j, x_s, y_s, z_s, x_u, y_u, z_u, x_n, y_n, z_n, x_mo, y_mo, z_mo, px_e, py_e, pz_e, px_v, py_v, pz_v, px_m, py_m, pz_m, px_ma, py_ma, pz_ma, px_j, py_j, pz_j, px_s, py_s, pz_s, px_u, py_u, pz_u, px_n, py_n, pz_n, px_mo, py_mo, pz_mo};

% Initial values

x_init_e = [-2.521092863852298E+10, 1.449279195712076E+11, -6.164888475164771E+05];
p_init_e = [-m_e*2.983983333368269E+04, -m_e*5.207633918704476E+03, m_e*6.169062303484907E-02];
x_init_v = [-1.075055502695123E+11, -3.366520720591562E+09, 6.159219802771119E+09];
p_init_v = [m_v*8.891598046362434E+02, -m_v*3.515920774124290E+04, -m_v*5.318594054684045E+02];
x_init_m = [-2.105262111032039E+10, -6.640663808353403E+10, -3.492446023382954E+09];
p_init_m = [m_m*3.665298706393840E+04, -m_m*1.228983810111077E+04, -m_m*4.368172898981951E+03];
x_init_ma = [2.079950549843747E+11, -3.143009598909217E+09, -5.178781274035379E+09];
p_init_ma = [m_ma*1.295003511847332E+03, m_ma*2.629442068920678E+04, m_ma*5.190097462980710E+02];
x_init_j = [5.989091025404842E+11, 4.391226571737888E+11, -1.523254483973965E+10];
p_init_j = [-m_j*7.901457175335077E+03, m_j*1.116228292421271E+04, m_j*1.306493838971390E+02];
x_init_s = [9.587061411378950E+11, 9.825654188805556E+11, -5.522074573564696E+10];
p_init_s = [-m_s*7.429888560624798E+03, m_s*6.738103977265501E+03, m_s*1.781267199919716E+02];
x_init_u = [2.158774693803164E+12, -2.054825232983900E+12, -3.562351187986338E+10];
p_init_u = [m_u*4.637499185676659E+03, m_u*4.627218902197726E+03, -m_u*4.289938933049831E+01];
x_init_n = [2.514853357481393E+12, -3.738847416615588E+12, 1.903951252367210E+10];
p_init_n = [m_n*4.466066630387768E+03, m_n*3.076399467848647E+03, -m_n*1.660621420555666E+02];
x_init_mo = [-2.552857888050620E+10, 1.446860363961675E+11, 3.593933517466486E+07];
p_init_mo = [-m_mo*2.927904627038706E+04, -m_mo*6.007566180814270E+03, -m_mo*1.577640655646029];

init_vals = [x_init_e, x_init_v, x_init_m, x_init_ma, x_init_j, x_init_s, x_init_u, x_init_n, x_init_mo, p_init_e, p_init_v, p_init_m, p_init_ma, p_init_j, p_init_s, p_init_u, p_init_n, p_init_mo]; 

h = 100000;

t_start = 0;

%t_end = 40000000;
%t_end = 139998000; % 2004 Venus transit
%t_end = 140100000; % after 2004 Venus transit
%t_end = 392261340; % 2012 Venus transit
%t_end = 392400000; % after 2012 Venus transit
%t_end = 3721949280; % 2117 Venus transit
%t_end = 3800000000; % after 2117 Venus transit
%t_end = 9974198460; % 2125 Venus transit
t_end = 11000000000; % after 2125 Venus transit

%t_end = 31558118; %Earth orbital time
%t_end = 19414166; %Venus orbital time
%t_end = 7600530; % Mercury orbital time

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

theta(end);
k = 1;

for i = 1:length(theta)
    if theta(i) < 5*10^-3
        index(k) = i;
        theta_small(k) = theta(i);
        k = k + 1;
    end
end

error_x = abs(Y(1,7) - Y(end,7));
error_y = abs(Y(1,8) - Y(end,8));
error_z = abs(Y(1,9) - Y(end,9));
feilvektor = [error_x error_y error_z]';

abs_feil = sqrt(error_x^2 + error_y^2 + error_z^2);


% Hamiltonian energy

H = 0;
for i = 1:length(Y(:,1))
    H(i) = hamiltonian_energy(Y(i,:),m);
end


% Plotting

h1 = animatedline; % earth line
 d = animatedline; % venus line
 m_1 = animatedline; % mercury line
 ma_1 = animatedline; % mars line
 q = animatedline; % jupiter line
 r = animatedline; % saturn line
 s = animatedline; % uranus line
 t1 = animatedline; % neptun line
 mo1 = animatedline; % moon line
 
 e = animatedline; % earth dot
 w = animatedline; % venus dot
 m_2 = animatedline; % mercury dot
 ma_2 = animatedline; % mars dot
 q2 = animatedline; % jupiter dot
 r2 = animatedline; % saturn dot
 s2 = animatedline; % uranus dot
 t2 = animatedline; % neptun dot
 mo2 = animatedline; % moon dot
 
 axis([-0.2*10^12,0.2*10^12,-0.2*10^12,0.2*10^12,-0.2*10^12,0.2*10^12])
 axis square
 
 hold on
 plot(0,0,'r*','Markersize',5)
 text(0,0,'sun')
 
 
 for k = 1:length(Y(:,1))
     addpoints(h1,Y(k,1),Y(k,2),Y(k,3));
     addpoints(d,Y(k,4),Y(k,5),Y(k,6));
     addpoints(m_1,Y(k,7),Y(k,8),Y(k,9));
     addpoints(ma_1,Y(k,10),Y(k,11),Y(k,12));
     addpoints(q,Y(k,13),Y(k,14),Y(k,15));
     addpoints(r,Y(k,16),Y(k,17),Y(k,18));
     addpoints(s,Y(k,19),Y(k,20),Y(k,21));
     addpoints(t1,Y(k,22),Y(k,23),Y(k,24));
     addpoints(mo1,Y(k,25),Y(k,26),Y(k,27));
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
     set(mo2,'Color','b','Marker','.','MarkerSize',10,'MaximumNumPoints',1);
     addpoints(mo2,Y(k,25),Y(k,26),Y(k,27));
     drawnow
 end
