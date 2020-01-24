function [ff, fx] = fill_ecdf(f, x)
    lb = min(x);
    ub = max(x);
    n = ub-lb + 1;
    fx = lb:ub;
    ff = zeros(1,n);
    pointer = 0;
    for i = 1:n
        if x(pointer+1)==fx(i)
            pointer = pointer + 1;
        end
        ff(i) = f(pointer);
    end

end



