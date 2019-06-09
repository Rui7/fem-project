function mat = subdomain(number)
%SUBDOMAIN Get material index
%   mat = subdomain(number)
%   where 1=Aluminium, 2=Steel, 3=Copper, 4=Electricity core
    if number == 1
        mat = 4;
    end
    if number == 2 || number == 4 || number == 6
        mat = 1;
    end
    if number == 3
        mat = 3;
    end
    if number == 5
        mat = 2;
    end
    