function y = BDF(y0, t0, dt, intTer, lambda, q, alpha)
%% This method use explicit method

    y = [y0];
    
    y_diff = [0];
    for time = 1:intTer

        f = source(t0+(time-1)*dt, alpha);
       
        if time == 1
            sol = (lambda+q(1))*y(end)+f;
            y_diff = [];
        else
            sol = (lambda+q(1))*y(end)-sum(y_diff.*q(2:time))+f;
        end
        y = [y, sol/q(1)];
        
        
        y_diff = [y(end)-y(end-1), y_diff];
    end
end

