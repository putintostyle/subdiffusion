function [alpha,result, max_result, sol] = subdiffu_main(alpha, order, Domain_size, Nx, dimension, T, pow_min, pow_max, modified)
    %% Domain_size
    
    Nx = Nx-2;
    dx = Domain_size/Nx;
    % dimension = 1;
    init = initial(Domain_size, dx, dimension);
    init = init(:);
    
    D = Laplacian(Nx, Nx, dx, dimension, 0);
    % T = 1;
    % dt = T/2^14;
    % alpha = 0.3;
    
    %%%%%
    % order = 2;
    % %%%%%
    ref = modif_subdiffusion(1, Nx, D, init, alpha, 1, 2^(-14),1);
    ref = ref(:, end);
    % 
    % x_grid = dx:dx:Domain_size;
    
    
    power = pow_min:pow_max;
    T_list = T./(2.^power);
    max_result = [];
    result = [];
    for time_s = length(T_list):-1:1
        if modified
            sol = modif_subdiffusion(order, Nx, D, init, alpha, T, T_list(time_s), 1);
        else
            sol = subdiffusion(order, Nx, D, init, alpha, T, T_list(time_s), 1);
        end
%         if time_s<length(T_list)
%             disp(size(sol_ref(:,1:2:end)))
%             disp(size(sol(:,1:end)))
%             if order == 1
%                 e = sum((sol-sol_ref(:,1:2:end)).^2, 1).^(1/2)*dx;
% %                 disp(e(1:6))
%             else
%                 e = sum((sol(:, 2:end)-sol_ref(:,1:2:end)).^2, 1).^(1/2)*dx;
%                 
%             end
%             max_result = [max_result, max(e)];
%         end
%         if order == 1
%             sol_ref = sol;
%         else
%             sol_ref = sol(:, 2:end);
%         end
        result = [result, sol(:,end)];
    end
    % surf(reshape(ref, [Nx, Nx]))
%     e = sum((result(:, 2:end)-result(:, 1:end-1)).^2, 1).^(1/2);
    e = sum((result-ref).^2, 1).^(1/2);
    disp(e)
    tmp = T_list(length(T_list):-1:1);
    A = [log(tmp(1:end)'), ones(length(e),1)];
    b = log(e(:));
    sol = A\b
end
