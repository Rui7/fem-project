function mat = subdomain(number)
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
    