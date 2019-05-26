function area=triangleArea(T,points)
A = points(:,T(1));
B = points(:,T(2));
C = points(:,T(3));
a = A-B;
b = A-C;
a = [a ; 0];
b = [b ; 0];
area = norm(cross(a,b))/2;

