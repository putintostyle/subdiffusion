function q = q_weight(time, dt, alpha)
    q = [];
    for i = time:-1:1
        q = [q, (time-i+1)^(1-alpha)-(time-i)^(1-alpha)];
    end
    q = q./(gamma(2-alpha)*(dt^alpha));
end



