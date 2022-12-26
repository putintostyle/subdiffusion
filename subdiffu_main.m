function [alpha,result, max_result, sol] = subdiffu_main(alpha, order, Domain_size, Nx, dimension, T, pow_min, pow_max)
    % Domain_size
    % = 1;
    
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
%     ref = subdiffusion(order, Nx, D, init, alpha, T, dt,1);
    % 
    % x_grid = dx:dx:Domain_size;
    % %%%%%
    % exact = 0;
    % for n=1:10
    %     exact = exact+16/pi/pi/pi*(ml(-n^2*pi^2*T^alpha, alpha, 1))/n/n/n*(1-(-1)^n)*sin(n*pi.*x_grid);
    % end
    % %%%%%
    
    power = pow_min:pow_max;
    T_list = T./(2.^power);
    max_result = [];
    result = [];
    for time_s = length(T_list):-1:1
        sol = subdiffusion(order, Nx, D, init, alpha, T, T_list(time_s), 1);
        
        if time_s<length(T_list)
            disp(size(sol_ref(:,1:2:end)))
            disp(size(sol))
            if order == 1
                e = sum(dx.*(sol-sol_ref(:,1:2:end)).^2, 1).^(1/2);
            else
                e = sum(dx.*(sol(:, 2:end)-sol_ref(:,1:2:end)).^2, 1).^(1/2);
            end
            max_result = [max_result, max(e)];
        end
        if order == 1
            sol_ref = sol;
        else
            sol_ref = sol(:, 2:end);
        end
        result = [result, sol(:,end)];
    end
    % surf(reshape(ref, [Nx, Nx]))
    e = sum((result(:, 2:end)-result(:, 1:end-1)).^2, 1).^(1/2);
    % e = sum((result-ref).^2, 1).^(1/2);
    tmp = T_list(length(T_list):-1:1);
    A = [log(tmp(2:end)'), ones(length(max_result),1)];
    b = log(max_result(:));
    sol = A\b
end
