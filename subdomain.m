function D = subdomain(number)
    k = [238, 20, 385, 1.6];
    if number == 1
        D = k(4);
    end
    if number == 2 || number == 4 || number == 6
        D = k(1);
    end
    if number == 3
        D = k(3);
    end
    if number == 5
        D = k(2);
    end
    