function [y, y_his, y_numerical] = parareal(y0, ts, te, dt, Dt, lambda, k, alpha, isParallel)
%     addpath('/Users/lvbingze/Desktop/subdissusion')
    coarseT = ts:Dt:te;
    
    M = length(coarseT);
    N = ceil(Dt/dt);
    
    qc = q_weight(M, Dt, alpha); %% [q(0), q(1...)]

    qf = q_weight(N, dt, alpha);
    
    qe = q_weight(ceil((te-ts)/dt)+1, dt, alpha);
%% record solutions in each iteration
    y_his = [];
%% main
    para_s = tic;
    yCoarse = BDF(y0, ts, Dt, M, lambda, qc, alpha);
    y_his = [y_his, yCoarse];
    yFine = yCoarse;
    
    for iteration = 1:k
        %% fine solver at each time-interval
        if isParallel
            parfor i = 2:M
                
                fineSol = BDF(yCoarse(i-1), coarseT(i-1), dt, N, lambda, qf, alpha);
                
                yFine(i) = fineSol(end);
            end
        else
            for i = 2:M
                
                fineSol = BDF(yCoarse(i-1), coarseT(i-1), dt, N, lambda, qf, alpha);
                
                yFine(i) = fineSol(end);
            end
        end
        %% coarse solver
        yCoarsediff = [];
        for i = 2:M

            f = source(ts+(i-2)*Dt, alpha);

            if i == 2 %% first step update y(t_1)
                sol = (qc(1)+lambda)*yCoarse(i-1)+f+[yFine(i)-yCoarse(i)];
                
            else
                
                sol = (qc(1)+lambda)*yCoarse(i-1)-sum(yCoarsediff.*qc(2:i-1))+f+[yFine(i)-yCoarse(i)];
                
            end
            
            yCoarse(i) = sol/qc(1);
            yCoarsediff = [yCoarse(i)-yCoarse(i-1), yCoarsediff];
            
        end
        y_his = [y_his; yCoarse];
    end
    
    y = yCoarse;
    disp(toc(para_s))
    %% Serial
    serial_s = tic;
    y_numerical = BDF(y0, ts, dt, ceil((te-ts)/dt)+1, lambda, qe, alpha);
    
    y_numerical = y_numerical(1:N:end);
    disp(toc(serial_s))
end