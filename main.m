Domain_size = 1;
Nx = 50-2;
dx = Domain_size/Nx;
dimension = 1;
init = initial(Domain_size, dx, dimension);
init = init(:);
D = Laplacian(Nx, Nx, dx, dimension, 1); %(nx, ny, h_step ,dim, method)
T = 1;
pow_2 = 12;
dt = 1/(2^pow_2);

alpha = 0.4;
eps = dx;
plot = false;
% [reference] = AllenCahn(1, Nx, D, init, alpha, T, dt, 1, 'l', eps); %(order, Nx, D, init, alpha, T, dt, dim, method, eps)

power = 10:15;
T_list = T./(2.^power);
my_data = [];
for alpha_list = [0.1, 0.4, 0.7, 0.9]
    result = [];
    for time_s = length(T_list):-1:1
        
        [sol] = AllenCahn(1, Nx, D, init, alpha, T, T_list(time_s), 1, 'l', eps);
        
        if time_s<length(T_list)
            disp(size(sol_ref(:,1:2:end)))
            disp(size(sol))
            e = sum((sol-sol_ref(:,1:2:end)).^2, 1).^(1/2);
            result = [result, max(e)];
        end
    sol_ref = sol;
    
    
    end
    my_data = [my_data; result];

end
save LTE my_data
% 
% e = sum((result-reference(:,end)).^2, 1).^(1/2);
tmp = T_list(length(T_list):-1:1)'; 
A = [log(tmp(1:end-1)), ones(length(result),1)];
b = log(result(:));
A\b





if plot
    if dimension == 2
    x = dx:dx:Domain_size;
    y = dx:dx:Domain_size;
    [xx, yy] = meshgrid(x, y);
    figure;
    subplot(2,2,2);
    surf(xx, yy, reshape(hist_arr(:,end), [Nx, Nx]));
    title('terminate')
    colorbar;
    caxis([-1 1]);
    axis equal;
    subplot(2,2,4);
    contourf(reshape(hist_arr(:,end), [Nx, Nx]));
    colorbar
    caxis([-1 1]);
    axis equal;
    
    subplot(2,2,1);
    surf(xx, yy, reshape(hist_arr(:,1), [Nx, Nx]));
    colorbar
    colorbar;
    caxis([-1 1]);
    title('initial')
    axis equal;
    subplot(2,2,3);
    contourf(reshape(hist_arr(:,1), [Nx, Nx]));
    colorbar
    colorbar;
    caxis([-1 1]);
    axis equal;
else
    x = dx:dx:Domain_size;
    figure;
    subplot(1,2,1);
    plot(x, reshape(hist_arr(:,1), [1,Nx]));
    title('initial')
   
    axis ([0 Domain_size -1 1]);
    subplot(1,2,2);
    plot(x, reshape(hist_arr(:,end), [1,Nx]));
    title('terminate')
    axis ([0 Domain_size -1 1]);
end

    
end
