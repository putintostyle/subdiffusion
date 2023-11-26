
function hist_arr = AllenCahn(method, Nx, D, init, alpha, T, dt, dim, eps)
   q = q_weight(T/dt, dt, alpha);
   % % % % % % % % % %         BDF1
   if method == 1 %linear
       hist_arr = [init(:)];
       if dim == 1
               M = sparse((q(1)+2).*eye(Nx)-eps^2.*D);
                    
       else
               M = sparse((q(1)+2).*eye(Nx*Nx)-eps^2.*D);
       end

       for iteration = progress(1:T/dt)
           if iteration == 1
               b = (q(1)+3).*hist_arr(:,end)-hist_arr(:,end).^3;
               root = M\b;
                
                    
           else
                diff_b = diff(hist_arr, 1, 2);

                q_r = q(iteration:-1:2);
                
                b = -1.*sum(diff_b.*q_r, 2)+(q(1)+3).*hist_arr(:,end)-hist_arr(:,end).^3;
                
                root  = M\b;
                
           end
           hist_arr = [hist_arr, root(:)];
       end
       
   elseif method == 2 % nonlinearfunction
        hist_arr = [init(:)];
        for iteration = progress(1:T/dt)
            options = optimset('Display','off');
            root = fsolve(@(U) polyU(U, hist_arr, q, D, iteration), hist_arr(:, end), options);
            hist_arr = [hist_arr, root(:)];
        end
        
    end
end


function F = polyU(U, hist_arr, q, D, iteration)
    % global q  D  iteration
    if iteration == 1
       sol = (q(1)+1).*hist_arr(:,end);
    elseif iteration > 1
       d_b = diff(hist_arr, 1, 2);

       d_q = q(iteration:-1:2);
       sol = -1.*sum(d_b.*d_q, 2)+(q(1)+1).*hist_arr(:,end);
    end
    % disp(size(D))
    % disp(size(U))
    diffu = eps^2.*D*U;
    for i = 1:length(U)
        F(i) = U(i)^3+q(1)*U(i)-diffu(i)-sol(i);
    end

end

