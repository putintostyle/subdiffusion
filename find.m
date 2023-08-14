function y = fine(y0, dt)
    yArr = BDF(y0, alpha, ts, te, dt, lambda);
    y = yArr(end);
end