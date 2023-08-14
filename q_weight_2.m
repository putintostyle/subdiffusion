function [q1,q2]  = q_weight_2(time, dt, alpha)
    g2 = gamma(2-alpha);
    g3 = gamma(3-alpha);
    q1 = [];
    q2 = [];
    for i = time:-1:1
        q1 = [q1, ((time-i+1)^(2-alpha)-(time-i)^(2-alpha))/g3-(time-i)^(1-alpha)/g2]; % q1(k) = q1^(n)_{n-k}
        q2 = [q2, ((time-i+1)^(2-alpha)-(time-i)^(2-alpha))/g3-(time-i+1)^(1-alpha)/g2];
    end
    q1 = q1./(dt^alpha);
    q2 = q2./(dt^alpha);
end



