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
thickness = 50;

n_elem = size(t,2);
n_nod = size(p,2);
elem_nod = 3;

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

dof = (1:n_nod)'; % rätt?? fattar inte dof

er = e([1 2 5],:);

conv_segments = [15 16];
conv_segments2 = [18 19];
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
end

coord = p';

q_n = @(T) alpha_c * (T - T_inf);
C_xy = @(ex,ey) [[1;1;1] ex ey];

B_bar = [0 1 0; 0 0 1];
%D= eye(2); %VAD ÄR D?

K = zeros(n_nod);
CC = zeros(n_nod);
f_l = zeros(n_nod,1);

for i=1:n_elem
    [ex,ey]=coordxtr(edof,coord,dof,elem_nod);
    C = C_xy(ex(i,:)',ey(i,:)');
    Ae = det(C)/2;
    D = k(subdomain(t(4,i)))*eye(2); %eller mer avancerad?
    Ke = C' \ B_bar' * D * B_bar / C * thickness * Ae;
    % Ke = flw2te(ex(i,:),ey(i,:),thickness,D);

    K = assem(edof(i,:),K,Ke); %finns snabbare assem i handledning
    
    Ce=plantml(ex(i,:),ey(i,:),1);
    CC=assem(edof(i,:),CC,Ce);
    
    fe = Q/3*thickness*Ae*[1;1;1];
    f_l = insert(edof(i,:),f_l,fe); %finns snabbare insert i handledning
end

K=sparse(K);
CC=sparse(CC);

%solve
% Tsnap=step1(K,C,d0,ip,f,pbound)

%DETTA HAR MICKE PRECIS LADDAT UPP
%Kefe
Ke = zeros(1,size(t,2));
fe = zeros(1,size(t,2));
for i = 1:size(t,2)
   triangle = t(:,i);
   ex = zeros(1,3);
   ey = zeros(1,3);
   eq = 0;
   for j = 1:3
    ex(j) = p(1,triangle(j));
    ey(j) = p(2,triangle(j));
   end
   D = k(subdomain(triangle(4)));
   if triangle(4) == 3
       eq = 100000;
   end
   %[Ke(i), fe(i)] = flw2te(ex, ey, thickness, D*eye(2), eq);
end
    
    