% load e, p and t
addpath(genpath('calfem/'))
load('e.mat')
load('p.mat')
load('t.mat')

% initialize
T0 = 25;
T_inf = 15;
Q = 1e5;
alpha_c = 100;
n_elem = size(t,2);
n_nod = size(p,2);
elem_nod = 3;
thickness = 50;

% material data
order = {'Aluminium'; 'Steel'; 'Copper'; 'Electricity core'};
E = [70, 210, 128, 500];
ny = [.33, .3, .36, .45];
alpha = [69e-6, 35e-6, 51e-6, 20e-6];
rho = [2710, 9700, 8930, 2000];
c_p = [903, 460, 386, 900];
k = [238, 20, 385, 1.6];

% t_areas = zeros(n_elem,1);
% for i = 1:n_elem
%     t_areas(i) = triangleArea(t(:,i),p);
% end

edof = (1:n_elem);
edof = [edof; t(1:3,:)]';

q_n = @(T) alpha_c * (T - T_inf);
C_xy = @(ex,ey) [[1;1;1] ex ey];

B_bar = [0 1 0; 0 0 1];
D= eye(2); %VAD ÄR D?

K = zeros(n_elem);
CC = zeros(n_elem);

for i=1:n_elem
    coord = p';
    dof = (1:n_nod)'; % rätt?? fattar inte dof
    [ex,ey]=coordxtr(edof,coord,dof,elem_nod);
    C = C_xy(ex(i,:)',ey(i,:)');
    Ae = det(C)/2;
    % D = k_material*eye(2) eller mer avancerad?
    Ke = C' \ B_bar' * D * B_bar / C * thickness * Ae;
    % Ke = flw2te(ex(i,:),ey(i,:),thickness,D);

    K = assem(edof(i,:),K,Ke); %finns snabbare assem i handledning
    
    Ce=plantml(ex(i,:),ey(i,:),1);
    CC=assem(edof(i,:),CC,Ce);
   
    
end

K=sparse(K);
C=sparse(C);

%solve
% Tsnap=step1(K,C,d0,ip,f,pbound)
