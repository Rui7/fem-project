function [ K_C,f_b ] = makeKC(K_C,f_b,edges_conv,coord, alpha_c, thickness, T_inf)
%MAKEKC Calculate the K_C matrix and the f_b vector
%   [ K_C,f_b ] = makeKC(K_C,f_b,edges_conv,coord, alpha_c, thickness, T_inf)

n_elem = size(edges_conv,2);

edgedof = 1:n_elem;
edgedof = [edgedof ; edges_conv]';

for i=1:n_elem
    p1=edges_conv(1,i);
    p2=edges_conv(2,i);
    y1=coord(p1,2);
    y2=coord(p2,2);
    L=(y2-y1);
    fe_ba = alpha_c*thickness*L/6 * [2 1;1 2];
    fe_b = alpha_c*T_inf*thickness*L*.5*[1;1];
    
    K_C = assem(edgedof(i,:),K_C,fe_ba);
    f_b = insert(edgedof(i,:),f_b,fe_b);
end
end

