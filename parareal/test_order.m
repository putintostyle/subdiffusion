y0 = 0;
ts = 0;
te = 1;
t = 6:12;
dt = 2.^(-t);
y_sol = [];
alpha = .5;
N = 2.^t+1;
lambda = -1;

for i = 1:length(dt)
    coarseT = 0:dt(i):1;
    exact = coarseT.^(3+alpha); 
    
    qe = q_weight(N(i), dt(i), alpha);
    
    y_numerical = BDF(y0, ts, dt(i), N(i), lambda, qe, alpha);
    length(exact)
    y_sol = [y_sol, max(abs(y_numerical-exact))];
end