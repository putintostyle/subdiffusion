function [order_list, error_list] =  test(alpha, order, Domain_size, Nx, dimension, T, pow_min, pow_max, ref_pow)
    Nx = Nx-2;
    dx = Domain_size/Nx;
    % dimension = 1;
    init = initial(Domain_size, dx, dimension);
    init = init(:);
    
    D = Laplacian(Nx, Nx, dx, dimension, 0);
    order_list = [];
    error_list = [];

    for a=alpha
        ref = subdiffusion(order, Nx, D, init, alpha, T, T/(2^ref_pow)/(10^4), 1);

        power = pow_min:pow_max;
        T_list = T./(2.^power)/(10^4);
        max_result = [];
        result = [];
        idx = 0;
        for time_s = length(T_list):-1:1

            sol = subdiffusion(order, Nx, D, init, alpha, T, T_list(time_s), 1);
            
            if time_s>0
                disp(size(ref(:,1:2^(ref_pow-power(end-idx)):end)))
                disp(size(sol(:,1:end)))
                if order == 1
                    e = sum((sol-ref(:,1:2^(ref_pow-power(end-idx)):end)).^2, 1).^(1/2)*dx;
    %                 disp(e(1:6))
                else
                    e = sum((sol(:, 2:end)-ref(:,1:2^(ref_pow-power(idx)):end)).^2, 1).^(1/2)*dx;
                    
                end
                max_result = [max_result, max(e)];
                idx = idx+1;
            end
%             if order == 1
%                 sol_ref = sol;
%             else
%                 sol_ref = sol(:, 2:end);
%             end
            result = [result, sol(:,end)];
        end
        tmp = T_list(length(T_list):-1:1);
        A = [log(tmp(1:end)'), ones(length(max_result),1)];
        b = log(max_result(:));
        sol = A\b
        order_list = [order_list, sol(1)];
        error_list = [error_list; max_result];
    end

end
