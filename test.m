function [order_list_1, order_list_2, error_list_1, error_list_2] =  test(alpha, order, Domain_size, Nx, dimension, T, pow_min, pow_max, ref_pow)
    Nx = Nx-2;
    dx = Domain_size/Nx;
    % dimension = 1;
    init = initial(Domain_size, dx, dimension);
    init = init(:);
    
    D = Laplacian(Nx, Nx, dx, dimension, 0);
    order_list_1 = [];
    order_list_2 = [];
    error_list_1 = [];
    error_list_2 = [];
    for a=alpha
        ref = subdiffusion(order, Nx, D, init, a, T, T/(2^ref_pow)./(10^4), 1);

        power = pow_min:pow_max;
        T_list = T./(2.^power)./(10^4);
        max_result_1 = [];
        max_result_2 = [];
        result = [];
        idx = 0;
        for time_s = length(T_list):-1:1

            sol = subdiffusion(order, Nx, D, init, a, T, T_list(time_s), 1);
            
%             if time_s>0
                disp(size(ref(:,1:2^(ref_pow-power(end-idx)):end)))
                disp(size(sol(:,1:end)))
                if order == 1
                    e_1 = sum((sol-ref(:,1:2^(ref_pow-power(end-idx)):end)).^2, 1).^(1/2)*dx;
                    if time_s<length(T_list)
                        e_2 = sum((sol-sol_ref(:,1:2:end)).^2, 1).^(1/2)*dx;
                    end
                    sol_ref = sol;
    %                 disp(e(1:6))
                else
                    solu = (sol(:, 2:end));
                    e_1 = sum((solu-ref(:,2:2^(ref_pow-power(end-idx)):end)).^2, 1).^(1/2)*dx;
                    if time_s<length(T_list)
                        e_2 = sum((solu-sol_ref(:,1:2:end)).^2, 1).^(1/2)*dx;
                    end
                    sol_ref = solu;
                end
                max_result_1 = [max_result_1, max(e_1)];
                if time_s<length(T_list)
                    max_result_2 = [max_result_2, max(e_2)];
                end
                idx = idx+1;
                
%             end
%             if order == 1
%                 sol_ref = sol;
%             else
%                 sol_ref = sol(:, 2:end);
%             end
            result = [result, sol(:,end)];
        end
        tmp = T_list(length(T_list):-1:1);
        A = [log(tmp(1:end)'), ones(length(max_result_1),1)];
        b = log(max_result_1(:));
        sol = A\b;
        order_list_1 = [order_list_1, sol(1)];
        error_list_1 = [error_list_1; max_result_1];
        A = [log(tmp(2:end)'), ones(length(max_result_2),1)];
        b = log(max_result_2(:));
        sol = A\b;
        order_list_2 = [order_list_2, sol(1)];
        error_list_2 = [error_list_2; max_result_2];
    end

end
