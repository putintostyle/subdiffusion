Domain_size = 1;
Nx = 100-2;
dx = Domain_size/Nx;
dimension = 1;
init = initial(Domain_size, dx, dimension);
init = init(:);
D = Laplacian(Nx, Nx, dx, dimension, 0);
T = 1;
dt = T/2^14;
alpha = 0.4;

%%%%%
order = 1;
% %%%%%
ref = subdiffusion(order, Nx, D, init, alpha, T, dt,1);
% 
% x_grid = dx:dx:Domain_size;
% %%%%%
% exact = 0;
% for n=1:10
%     exact = exact+16/pi/pi/pi*(ml(-n^2*pi^2*T^alpha, alpha, 1))/n/n/n*(1-(-1)^n)*sin(n*pi.*x_grid);
% end
% %%%%%

power = 7:12;
T_list = T./(2.^power);
result = [];
for time_s = 1:length(T_list)
    sol = subdiffusion(order, Nx, D, init, alpha, T, T_list(time_s), 1);
    result = [result, sol];
end
% surf(reshape(ref, [Nx, Nx]))
% e = sum((result(:, 2:end)-result(:, 1:end-1)).^2, 1).^(1/2);
e = sum((result-ref).^2, 1).^(1/2);
tmp = T_list(:); 
A = [log(tmp(:)), ones(length(e),1)];
b = log(e(:));
A\b  
