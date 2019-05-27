% load e, p and t
clear all;
addpath(genpath('calfem/'))
load('e.mat')
load('p.mat')
load('t.mat')

p=p*1e-3;

% initialize
T0 = 25;
T_inf = 15;
Q = 1e5;
alpha_c = 100;
thickness = 50e-3;

n_elem = size(t,2);
n_nod = size(p,2);
elem_nod = 3;

% material data
order = {'Aluminium'; 'Steel'; 'Copper'; 'Electricity core'};
E = [70, 210, 128, 500]*1e9;
ny = [.33, .3, .36, .45];
alpha = [69e-6, 35e-6, 51e-6, 20e-6];
rho = [2710, 9700, 8930, 2000];
c_p = [903, 460, 386, 900];
k = [238, 20, 385, 1.6];

% t_areas = zeros(n_elem,1);
% for i = 1:n_elem
%     t_areas(i) = triangleArea(t(:,i),p);
% end

a0 = T0*ones(n_nod,1);

edof = (1:n_elem);
edof = [edof; t(1:3,:)]';

dof = (1:n_nod)'; % rätt?? fattar inte dof

er = e([1 2 5],:);

coord = p';

q_n = @(T) alpha_c * (T - T_inf);
C_xy = @(ex,ey) [[1;1;1] ex ey];

B_bar = [0 1 0; 0 0 1];

K = zeros(n_nod);
CC = zeros(n_nod);
f_l = zeros(n_nod,1);

[ex,ey]=coordxtr(edof,coord,dof,elem_nod);

for i=1:n_elem
    mat_index = subdomain(t(4,i)); % index of material constants
    rhoe = rho(mat_index);
    ce = c_p(mat_index);
    D = k(mat_index)*eye(2); %eller mer avancerad?
    
    eq = 0;
    C = C_xy(ex(i,:)',ey(i,:)');
    Ae = det(C)/2;%*1e-6;
    
    if mat_index == 3
       eq = Q/3*thickness*Ae;
    end
    
    %Ke = C' \ B_bar' * D * B_bar / C * thickness * Ae;
    [Ke,fe] = flw2te(ex(i,:),ey(i,:),thickness,D,eq);

    K = assem(edof(i,:),K,Ke); %finns snabbare assem i handledning
    
    Ce=plantml(ex(i,:),ey(i,:),rhoe*ce);
    CC=assem(edof(i,:),CC,Ce);
    
    %fe = Q/3*thickness*Ae*[1;1;1];
    f_l = insert(edof(i,:),f_l,fe); %finns snabbare insert i handledning
end

K=sparse(K);
CC=sparse(CC);

conv_segments_al = [15,18];
conv_segments_st = [16,19];
%conv_segments2 = [18 19]; 
edges_conv_al = [];
edges_conv_st = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments_al)
        edges_conv_al = [edges_conv_al er(1:2,i)];
    elseif ismember(er(3,i),conv_segments_st)
        edges_conv_st = [edges_conv_st er(1:2,i)];
    end
end

K_C = zeros(n_nod);
f_b = zeros(n_nod,1);

[K_C,f_b] = makeKC(K_C,f_b,edges_conv_al,coord, alpha_c, thickness, T_inf);
[K_C,f_b] = makeKC(K_C,f_b,edges_conv_st,coord, alpha_c, thickness, T_inf);

Kprim = sparse(K+K_C);
f = sparse(f_l + f_b);

% bc = [edges_conv_al(:); edges_conv_st(:)];
% bc = unique(bc);
% bc = [bc T_inf*ones(length(bc),1)];

%solve
a=solveq(Kprim,f);
% Tsnap=step1(K,CC,a0,ip,f,pbound)

%plot
eT=extract(edof,a);
figure()
fill(ex',ey',eT','EdgeColor','none')
title('Temperature distribution [C]')
colormap(hot);
colorbar;
xlabel('x-position [mm]')
ylabel('y-position [mm]')
%axis equal

%DETTA HAR MICKE PRECIS LADDAT UPP
%Kefe
% Ke = zeros(1,size(t,2));
% fe = zeros(1,size(t,2));
% for i = 1:size(t,2)
%    triangle = t(:,i);
%    ex = zeros(1,3);
%    ey = zeros(1,3);
%    eq = 0;
%    for j = 1:3
%     ex(j) = p(1,triangle(j));
%     ey(j) = p(2,triangle(j));
%    end
%    D = k(subdomain(triangle(4)));
%    if triangle(4) == 3
%        eq = 100000;
%    end
%    %[Ke(i), fe(i)] = flw2te(ex, ey, thickness, D*eye(2), eq);
% end
    
% edgedof = 1:size(er,2);
% edgedof = [edgedof ; er(1:2,:)];
    