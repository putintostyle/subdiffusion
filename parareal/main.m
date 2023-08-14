clear all
y0 = 0;
ts = 0;
te = 1;
dt = 1e-5;
Dt = dt*10;
lambda = -1;
k = 10;
alpha = .5;
coarseT = ts:Dt:te;
exact = coarseT.^(3); 
isParallel = 1;

 
[y, y_his, y_num] = parareal(y0, ts, te, dt, Dt, lambda, k, alpha, isParallel);

%%
res = [];
parfor i=2:k
    res(i-1) = max(abs(y_his(i,:)-exact));
end
close all
figure
semilogy(2:k, (res), 'b')
hold on
scatter(2:k, (res), '^')

xlabel('iteration')
ylabel('max absolue error to exact solution')
hold off
%%
close all

figure
for i=1:k
loglog(coarseT, abs(y_his(i, :)-exact), '-')
hold on
end
% y_1 = abs(y_num-exact);
loglog(coarseT, abs(y_num-exact), 'r')
xlabel('time')
ylabel('absolue error to exact solution')
% loglog(coarseT, exact, 'r')
hold off


