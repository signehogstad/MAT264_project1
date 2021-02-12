hold off
clear
clc

% List of bodies being considered
% 1) Earth
% 2) Venus
% 3) Mercury
% 4) Mars
% 5) Jupiter
% 6) Saturn
% 7) Uranus
% 8) Neptun
% 9) Moon
% 10) Sun

% Constants
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
m_su = 1.98847*10^30;

m = [m_e m_v m_m m_ma m_j m_s m_u m_n m_mo m_su];


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

x_su = @(px_su) position(px_su,m_su);
y_su = @(py_su) position(py_su,m_su);
z_su = @(pz_su) position(pz_su,m_su);


% Momentum vectors

px_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_e = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Earth',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_v = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Venus',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_m = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Mrcry',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_ma = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Marss',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_j = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Jupit',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_s = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Satur',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_u = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Uranu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_n = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Neptu',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_mo = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Mooon',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_mo = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Mooon',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_mo = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Mooon',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

px_su = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Suuun',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,1);
py_su = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Suuun',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,2);
pz_su = @(x_e,y_e,z_e,x_v,y_v,z_v,x_m,y_m,z_m,x_ma,y_ma,z_ma,x_j,y_j,z_j,x_s,y_s,z_s,x_u,y_u,z_u,x_n,y_n,z_n,x_mo,y_mo,z_mo,x_su,y_su,z_su) momentum_ssb('Suuun',[x_e y_e z_e x_v y_v z_v x_m y_m z_m x_ma y_ma z_ma x_j y_j z_j x_s y_s z_s x_u y_u z_u x_n y_n z_n x_mo y_mo z_mo x_su y_su z_su],m,3);

funcs = {x_e, y_e, z_e, x_v, y_v, z_v, x_m, y_m, z_m, x_ma, y_ma, z_ma, x_j, y_j, z_j, x_s, y_s, z_s, x_u, y_u, z_u, x_n, y_n, z_n, x_mo, y_mo, z_mo, x_su, y_su, z_su, px_e, py_e, pz_e, px_v, py_v, pz_v, px_m, py_m, pz_m, px_ma, py_ma, pz_ma, px_j, py_j, pz_j, px_s, py_s, pz_s, px_u, py_u, pz_u, px_n, py_n, pz_n, px_mo, py_mo, pz_mo, px_su, py_su, pz_su};

% Initial values

x_init_e = [-2.627892928682480E+10, 1.445102393586391E+11, 3.022818135935813E+07];
p_init_e = [-m_e*2.983052803283506E+04, -m_e*5.220465685407924E+03, -m_e*1.014621798034465E-01];
x_init_v = [-1.085735509178141E+11, -3.784200933160055E+09, 6.190064472977990E+09];
p_init_v = [m_v*8.984651054838754E+02, -m_v*3.517203950794635E+04, -m_v*5.320225582712421E+02];
x_init_m = [-2.212062175862221E+10, -6.682431829610253E+10, -3.461601353176080E+09];
p_init_m = [m_m*3.666229236478603E+04, -m_m*1.230266986781422E+04, -m_m*4.368336051784789E+03];
x_init_ma = [2.069270543360729E+11, -3.560689811477710E+09, -5.147936603828508E+09];
p_init_ma = [m_ma*1.304308812694964E+03, m_ma*2.628158892250333E+04, m_ma*5.188465934952333E+02];
x_init_j = [5.978411018921824E+11, 4.387049769612203E+11, -1.520170016953275E+10];
p_init_j = [-m_j*7.892151874487445E+03, m_j*1.114945115750926E+04, m_j*1.304862310943014E+02];
x_init_s = [9.576381404895931E+11, 9.821477386679871E+11, -5.518990106544006E+10];
p_init_s = [-m_s*7.420583259777167E+03, m_s*6.725272210562052E+03, m_s*1.779635671891335E+02];
x_init_u = [2.157706693154863E+12, -2.055242913196469E+12, -3.559266720965648E+10];
p_init_u = [m_u*4.646804486524291E+03, m_u*4.614387135494279E+03, -m_u*4.306254213333638E+01];
x_init_n = [2.513785356833091E+12, -3.739265096828156E+12, 1.907035719387865E+10];
p_init_n = [m_n*4.475371931235400E+03, m_n*3.063567701145199E+03, -m_n*1.662252948584049E+02];
x_init_mo = [-2.659657952880802E+10, 1.442683561835989E+11, 6.678400538153946E+07];
p_init_mo = [-m_mo*2.926974096953943E+04, -m_mo*6.020397947517719E+03, -m_mo*1.740793458484102];
x_init_su = [-1.068000648301820E+09, -4.176802125684930E+08, 3.084467020687090E+07];
p_init_su = [m_su*9.305300847631915E+00, -m_su*1.283176670344807E+01, -m_su*1.631528028381386E-01];

init_vals = [x_init_e, x_init_v, x_init_m, x_init_ma, x_init_j, x_init_s, x_init_u, x_init_n, x_init_mo, x_init_su, p_init_e, p_init_v, p_init_m, p_init_ma, p_init_j, p_init_s, p_init_u, p_init_n, p_init_mo, p_init_su]; 

h = 5000;

t_start = 0;

%t_end = 40000000;
%t_end = 139998000; % 2004 Venus transit
t_end = 150100000; % after 2004 Venus transit
%t_end = 392261340; % 2012 Venus transit
%t_end = 392400000; % after 2012 Venus transit
%t_end = 3721949280; % 2117 Venus transit
%t_end = 3800000000; % after 2117 Venus transit
%t_end = 3974198460; % 2125 Venus transit
%t_end = 4000000000; % after 2125 Venus transit
%t_end = 7808527980; % 2247 Venus transit
%t_end = 8060791080; % 2255 Venus transit
%t_end = 8100000000; % after 2255 Venus transit

%t_end = 31558118; %Earth orbital time
%t_end = 19414166; %Venus orbital time
%t_end = 7600530; % Mercury orbital time

% Try to solve

[Y,t] = euler_forward_ssb(funcs, init_vals, h, t_start, t_end);


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
 su1 = animatedline; % sun line
 
 e = animatedline; % earth dot
 w = animatedline; % venus dot
 m_2 = animatedline; % mercury dot
 ma_2 = animatedline; % mars dot
 q2 = animatedline; % jupiter dot
 r2 = animatedline; % saturn dot
 s2 = animatedline; % uranus dot
 t2 = animatedline; % neptun dot
 mo2 = animatedline; % moon dot
 su2 = animatedline; % sun dot
 
 axis([-0.2*10^12,0.2*10^12,-0.2*10^12,0.2*10^12,-0.2*10^12,0.2*10^12])
 axis square
 
 hold on
 
 
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
     addpoints(su1,Y(k,28),Y(k,29),Y(k,30));
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
     set(su2,'Color','y','Marker','.','MarkerSize',30,'MaximumNumPoints',1);
     addpoints(su2,Y(k,28),Y(k,29),Y(k,30));
     drawnow
 end