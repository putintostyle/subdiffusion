function [order_list, error_list] = batch(alpha, order, Domain_size, Nx, dimension, T, pow_min, pow_max)
    order_list = [];
    error_list = [];

    for a=alpha
        [alpha_a,result, max_result, sol] = subdiffu_main(a, order, Domain_size, Nx, dimension, T, pow_min, pow_max);
        error_list = [error_list; result];
        order_list = [order_list, sol(1)];
    end

end
