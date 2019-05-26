% load e, p and t
load('e.mat')
load('p.mat')
load('t.mat')

T0 = 25;
T_inf = 15;
Q = 1e5;
a_c = 100;

order = {'Aluminium'; 'Steel'; 'Copper'; 'Electricity core'};
E = [70, 210, 128 500];
ny = [.33, .3, .36, .45];
alpha = [69e-6, 35e-6, 51e-6, 20e-6];
rho = [2710, 9700, 8930, 2000];
c_p = [903, 460, 386, 900];
k = [238, 20, 385, 1.6];



