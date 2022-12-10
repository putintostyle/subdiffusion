function init = initial(Domain_size, dx, dim)
    if dim ==1
        x = dx:dx:Domain_size;
        init = 0.05.*sin(x);
    else
        x = dx:dx:Domain_size;
        y = dx:dx:Domain_size;
        [xx, yy] = meshgrid(x, y);
        init = sin(xx).*sin(yy);
%     surf(u)
    end
end
