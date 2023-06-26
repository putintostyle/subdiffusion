
function y = parareal(y0, ts, te, dt, Dt, lambda, k, alpha)
    
    coarseT = ts:Dt:te;
    
    M = length(coarseT);
    N = ceil(Dt/dt);

    qc = q_weight(M, Dt, alpha); %% [q(0), q(1...)]

    qf = q_weight(N, dt, alpha);

%% record solutions in each iteration
    y_his = [];
%% main

    yCoarse = BDF(y0, ts, Dt, M, lambda, qc, alpha);
    y_his = [y_his, yCoarse];
    yFine = yCoarse;
    
    for iteration = 1:k
        for i = 2:M
            fineSol = BDF(yCoarse(i-1), coarseT(i-1), dt, N, lambda, qf, alpha);
            
            yFine(i) = fineSol(end);
        end
        
        yCoarsediff = [0];
        for i = 2:M

            f = source(ts+(i-2)*Dt, alpha);

            if i == 2 %% first step update y(t_1)
                sol = (qc(1)+lambda)*yCoarse(i-1)+f+yFine(i)-yCoarse(i);
                yCoarsediff = [];
            else
                yCoarsediff
                sol = (qc(1)+lambda)*yCoarse(i-1)-sum(yCoarsediff.*qc(2:i-1))+f+yFine(i)-yCoarse(i);
                
            end
            
            yCoarse(i) = sol/qc(1);
            yCoarsediff = [yCoarse(i)-yCoarse(i-1), yCoarsediff];
            
        end
    end
    y = yCoarse;
end