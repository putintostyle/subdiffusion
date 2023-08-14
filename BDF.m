function y = BDF(y0, alpha, ts, te, dt, lambda)
    y = [y0];
    intTer = ceil((te-ts)/dt);
    q = q_weight(intTer, dt, alpha);
    ydiff = [];
    for time = 1:(intTer+1)
        if time == 1
            y = [y, y0*q(1)/(q(1)-lambda)];
            ydiff = [ydiff, y(2)-y0];
        else
            sol = (1*y(end-1)*q(1)-ydiff.*q(2:time))/(q(1)-lambda);
            y = [y, sol];
            ydiff = [ydiff, y(end)-y(end-1)];
        end
    end
end

