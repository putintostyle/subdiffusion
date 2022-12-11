function [hist_arr] = AllenCahn(order, Nx, D, init, alpha, T, dt, dim, method, eps)
   q = q_weight(T/dt, dt, alpha);
   % % % % % % % % % %         BDF1
   if order == 1
       hist_arr = [init(:)];
       for iteration = progress(1:T/dt)
            if method == 'c'
                fun = @(U)polyU_1(U, hist_arr, q, D, iteration);
                options = optimoptions('fsolve','Algorithm','levenberg-marquardt');

                x_0 = hist_arr(:,end);

                
                root = fsolve(fun, x_0, options);
            elseif method == 'l'
                b = (q(1)+3).*hist_arr(:,end)-hist_arr(:,end).^3;
            
            
                if iteration == 1
                    if dim == 1
                        root  = linsolve(((q(1)+2).*eye(Nx)-eps^2.*D),(b));
                    elseif dim ==2
                        root  = linsolve(((q(1)+2).*eye(Nx^2)-eps^2.*D),(b));
                    end
                    
                    
                else
                    for time = 1:iteration-1
                         
                        b = b-q(time+1).*(hist_arr(:,end-(time-1))-hist_arr(:,end-(time-1)-1));     
                             
                        if dim == 1
                            root  = linsolve(((q(1)+2).*eye(Nx)-eps^2.*D),(b));
                        elseif dim ==2
                            root  = linsolve(((q(1)+2).*eye(Nx^2)-eps^2.*D),(b));
                        end     
                    end
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

