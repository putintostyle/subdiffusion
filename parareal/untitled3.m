clear all

y0 = 0;
ts = 0;
te = 1;
dt = 1e-5;
Dt = dt*100;
lambda = -1;
k = 10;
alpha = .5;
coarseT = ts:Dt:te;
exact = coarseT.^(3); 
N = ceil((te-ts)/dt)+1;
qe = q_weight(N, dt, alpha);
y_numerical = BDF(y0, ts, dt, N, lambda, qe, alpha);

%%
close all
figure
% for i=1:k
% semilogx(coarseT, y_his(i, :), 'o')
% hold on
% end
hold on
semilogx(0:dt:te, y_numerical, 'bx')
semilogx(coarseT, exact, 'r')
hold off


