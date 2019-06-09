function f0 = makef0(ex,ey,D,nu,alpha,t,T,T0)
%UNTITLED2 Summary of this function goes here
%   f0 = makef0(ex,ey,D,nu,alpha,t,T,T0)

    Ae=1/2*det([ones(3,1) ex' ey']);
    
    C=[ 1  ex(1) ey(1)   0          0          0  
        0         0        0   1   ex(1)   ey(1)
        1  ex(2) ey(2)   0          0          0  
        0         0        0   1   ex(2)   ey(2)
        1  ex(3) ey(3)   0          0          0  
        0         0        0   1   ex(3)   ey(3)];

    B=[0 1 0 0 0 0
       0 0 0 0 0 1
%      0 0 0 0 0 0
       0 0 1 0 1 0]*inv(C);
   
   colD=size(D,2);
       if colD>3
         Cm=inv(D);
         Dm=inv(Cm([1 2 4],[1 2 4]));
       else
         Dm=D;
       end

    f0=B'*Dm*(1+nu)*alpha*t*1/3*([1,1,1]*(T-T0))*Ae*[1,1,0]';

end

