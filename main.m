% load e, p and t
clear all;
addpath(genpath('calfem/'))
load('e.mat')
load('p.mat')
load('t.mat')

plot=1; %ändra till 1 om du vill plotta

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
    
    Ce=plantml(ex(i,:),ey(i,:),rhoe*ce*thickness);
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
f=Kprim*a;

dt=3600*2/4;
T=3600*10/4;
alpha_method=1; %implicit
nsnap = 4;
nhist=1;
time=3600*[2,4,6,8]/4;

ip=[dt,T,alpha_method,[nsnap, nhist, time, dof']];
pbound=[];

Tsnap=step1(Kprim,CC,a0,ip,f,pbound);

eT=extract(edof,a);
%maxT=max(max(full(eT)))

maxT=zeros(nsnap,1);
% 
% step_size=.5*3600;
% init=.5*3600;
% steps=3;
% time=init:step_size:init+steps*step_size;
% 
% 
% next_step=step_size*Kprim+CC;
% prev_step=CC*a+step_size*f;
% 
% an=solveq(next_step,prev_step)

if plot == 1 % ändra högst upp om du vill plotta!
    %plot
    figure()
    fill(ex',ey',eT','EdgeColor','none')
    title('Temperature distribution [C]')
    colormap(hot);
    colorbar;
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    %axis equal
    
    Tmin=min(min(Tsnap));
    Tmax=max(max(Tsnap));
    figure()
    for i=1:nsnap
        eT=extract(edof,Tsnap(:,i));
        maxT(i) = max(max(eT));
        subplot(2,2,i)
        fill(ex',ey',eT','EdgeColor','none')
        title(['Temperature distribution [C] at T=',num2str(time(i)/3600),'h'])
        colormap(hot);
        colorbar;
        xlabel('x-position [m]')
        ylabel('y-position [m]')
        axis([0 .025 0 .05])
        %caxis([Tmin Tmax])
    end
end

%von Mises

K = zeros(2*n_nod);
f0 = zeros(2*n_nod,1);

ptype = 2; %plane strain
ep = [ptype thickness];

edof_dis = (1:n_elem);
tt=[t(1:3,:)*2-1; t(1:3,:)*2];
edof_dis = [edof_dis; tt]';
%edof_dis2
edof_dis = zeros(n_elem,7);
edof_dis(:,1) = (1:n_elem);
for i=1:n_elem
    edof_dis(i,2:end) = [t(1,i)*2-1, t(1,i)*2, t(2,i)*2-1, t(2,i)*2, t(3,i)*2-1, t(3,i)*2];
end
edof_dis

for i=1:n_elem
    
    mat_index = subdomain(t(4,i)); % index of material constants
    D = hooke(ptype,E(mat_index),ny(mat_index));
    nu=ny(mat_index);
    alphae =alpha(mat_index);
    
    Ke=plante(ex(i,:),ey(i,:),ep,D);
    K = assem(edof_dis(i,:),K,Ke); %finns snabbare assem i handledning
    
    nodes=t(1:elem_nod,i);
    T=full(a(nodes));
    
    fe0=makef0(ex(i,:),ey(i,:),D,nu,alphae,thickness,T,T0);
    f0 = insert(edof_dis(i,:),f0,fe0); %finns snabbare insert i handledning
end

K=sparse(K);
f0=sparse(f0);

conv_segments_ux = [14,17 ];
conv_segments_uy = [8,9,12,13]; 
bc = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments_ux)
        bc_el=[er(1:2,i)*2-1 [0;0]];
        bc=[bc;bc_el];
    elseif ismember(er(3,i),conv_segments_uy)
        bc_el=[er(1:2,i)*2 [0;0]];
        bc=[bc;bc_el];
    end
end
bc=unique(bc,'rows');

%få ut a-vektorn när vi räknat ut f0
a_dis=solveq(K,f0,bc);
ed = extract(edof_dis,a_dis);
   
%Seff=zeros();
s=zeros(n_elem,4);
for i=1:n_elem
    mat_index = subdomain(t(4,i)); % index of material constants
    
    %beräkna sigma och epsilon
    [es,et]=plants(ex(i,:),ey(i,:),ep,D,ed);
    s = s+es;
    
    %beräkna von mises
    h=0;
    %st=[0 0 0];
    %est=0;
    %mp=[E(mat_index),ny(mat_index),h];
    %[es,deps,st] = mises(ptype,mp,est,st);
end

vonmises = @(es) sqrt(es(:,1).*es(:,1) + es(:,2).*es(:,2) + es(:,3).*es(:,3) ... 
    - es(:,1).*es(:,2) -  es(:,1).*es(:,3) - es(:,2).*es(:,3) + 3* es(:,4).*es(:,4));

eff=vonmises(s);
eff_nod = zeros(size(p,2),1);
for i = 1:size(eff_nod)
    [c0 , c1] = find(edof(:,2:4)==i);
    eff_nod(i,1) = sum(eff(c0)/size(c0,1));
end

eff_e = extract(edof, eff_nod);

figure
fill(ex',ey',eff_e')%,'EdgeColor','none')
colorbar
title(['Von Mises Stresses at Stationary with T_{\infty}= ',num2str(T_inf), '°C'])
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis([0 .025 0 .05])
[Y,I] = max(eff);
I

max(max(eff_e))
min(min(eff_e))

%plotting displacements
pdis = zeros(size(p));
udisx = a_dis(1:2:end);
udisy = a_dis(2:2:end);
<<<<<<< HEAD
magnitude = 100;
=======
magnitude = 1000;
>>>>>>> 0069d95dd8bca45d55698b5a1c6321e23905f8a6
pdis(1,:) = p(1,:)+udisx'*magnitude;
pdis(2,:) = p(2,:)+udisy'*magnitude;
coorddis = pdis';
[exdis,eydis]=coordxtr(edof,coorddis,dof,elem_nod);
figure
title(['Displacement field with T_{\infty}= ',num2str(T_inf), '°C'])
fill(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
hold on
fill(exdis',eydis',[0 0 0],'EdgeColor','none',   'FaceAlpha', 0.3);
title(['Displacement field with T_{\infty}= ',num2str(T_inf), '°C'])
xlabel('x-position [m]')
ylabel('y-position [m]')
margin = 0.01;
axis([0-margin .025+margin 0-margin .05+margin])