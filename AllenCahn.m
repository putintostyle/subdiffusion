function hist_arr = AllenCahn(order, Nx, D, init, alpha, T, dt, dim, eps)
   q = q_weight(T/dt, dt, alpha);
   % % % % % % % % % %         BDF1
   if order == 1
       hist_arr = [init(:)];
       if dim == 1
               M = sparse((q(1)+2).*eye(Nx)-eps^2.*D);
                    
       else
               M = sparse((q(1)+2).*eye(Nx*Nx)-eps^2.*D);
       end
       for iteration = progress(1:T/dt)
           if iteration == 1
               b = (q(1)+3).*hist_arr(:,end)-hist_arr(:,end).^3;
               root  = M\b;
          
                    
                    
           else
                diff_b = diff(hist_arr, 1, 2);

                q_r = q(iteration:-1:2);
                
                b = -1.*sum(diff_b.*q_r, 2)+(q(1)+3).*hist_arr(:,end)-hist_arr(:,end).^3;
                
                root  = M\b;
                
           end
           hist_arr = [hist_arr, root(:)];
       end
       
   elseif oder == 2
       history_arr = [zeros(Nx,1),zeros(Nx,1)];
   end

% function F = polyU_1(U, hist_arr, q, D, iteration)
% 
%     for i = 1:length(U)
%         F(i) = U(i)^3+q(1)*U(i);
%     end
%     b = 1.*hist_arr(:,end);
%     if iteration > 1
%        for time = 1:iteration-1
%                          
%            b = b-q(time+1).*(hist_arr(:,end-(time-1))-hist_arr(:,end-(time-1)-1));     
%        end
% 
%     else
%         b = b+hist_arr(:,end).*q(1);
%     end
% 
%    F = F+eps^2.*D-b;

