function [hist_arr] = AllenCahn(order, Nx, D, init, alpha, T, dt, dim, method, eps)
   q = q_weight(T/dt, dt, alpha);
   % % % % % % % % % %         BDF1
   if order == 1
       hist_arr = [init(:)];
       if dim == 1
               M = sparse((-q(1).*eye(Nx)+D));
                    
       else
               M = sparse((-q(1).*eye(Nx*Nx)+D));
       end
       for iteration = progress(1:T/dt)
            if method == 'c'
                fun = @(U)polyU_1(U, hist_arr, q, D, iteration);
                options = optimoptions('fsolve','Algorithm','levenberg-marquardt');

                x_0 = hist_arr(:,end);

                
                root = fsolve(fun, x_0, options);
            elseif method == 'l'

                b = (q(1)+3).*hist_arr(:,end)-hist_arr(:,end).^3;
            
            
                if iteration == 1
                  
                    root  = M\b;
                end
                    
                    
                else
                    for time = 1:iteration-1
                         
                        b = b-q(time+1).*(hist_arr(:,end-(time-1))-hist_arr(:,end-(time-1)-1));     
                             
                        
                    
                    end
                    root  = M\b;
                end
            end
            hist_arr = [hist_arr, root(:)];
       end
   elseif oder == 2
       history_arr = [zeros(Nx,1),zeros(Nx,1)];
   end

function F = polyU_1(U, hist_arr, q, D, iteration)

    for i = 1:length(U)
        F(i) = U(i)^3+q(1)*U(i);
    end
    b = 1.*hist_arr(:,end);
    if iteration > 1
       for time = 1:iteration-1
                         
           b = b-q(time+1).*(hist_arr(:,end-(time-1))-hist_arr(:,end-(time-1)-1));     
       end

    else
        b = b+hist_arr(:,end).*q(1);
    end

   F = F+eps^2.*D-b;

