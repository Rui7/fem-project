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
T_inf = 25; %uppgift b
Q = 1e5;
%Q = 1e5*1.6^2; %uppgift a, del2
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

a0 = T0*ones(n_nod,1);

edof = (1:n_elem);
edof = [edof; t(1:3,:)]';

dof = (1:n_nod)';

er = e([1 2 5],:);

coord = p';

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
    D = k(mat_index)*eye(2);
    
    eq = 0;
    C = C_xy(ex(i,:)',ey(i,:)');
    Ae = det(C)/2;
    
    if mat_index == 4
       eq = Q;
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
dt=3600*2/4;
T=3600*8/4;
alpha_method=1; %implicit
nsnap = 4;
nhist=1;
time=3600*[0,2,4,6]/4;

ip=[dt,T,alpha_method,[nsnap, nhist, time, dof']];
pbound=[];

Tsnap=step1(Kprim,CC,a0,ip,f,pbound);

Tmin=min(min(Tsnap));
Tmax=max(max(Tsnap));
figure()
for i=1:nsnap
    eT=extract(edof,Tsnap(:,i));
    subplot(2,2,i)
    fill(ex',ey',eT','EdgeColor','none')
    title('Temperature distribution [C]')
    colormap(hot);
    colorbar;
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    axis([0 .025 0 .05])
    caxis([Tmin Tmax])
end

%plot
eT=extract(edof,a);
figure()
fill(ex',ey',eT','EdgeColor','none')
title('Temperature distribution [C]')
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
%axis equal

    