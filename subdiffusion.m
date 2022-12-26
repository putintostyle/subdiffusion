function history_arr = subdiffusion(order, Nx, D, init, alpha, T, dt, dim)
    
    q = q_weight(T/dt, dt, alpha);
    if dim == 1
        if order == 1
            history_arr = [zeros(Nx,1)];
        else
            history_arr = [zeros(Nx,1),zeros(Nx,1)];
        end
    else
        history_arr = [zeros(Nx*Nx,1),zeros(Nx*Nx,1)];
    end
    init = init(:);
    source = D*init;
    if order == 1
       if dim == 1
               M = sparse((-q(1).*eye(Nx)+D));
                    
       else
               M = sparse((-q(1).*eye(Nx*Nx)+D));
       end

    elseif order == 2
       if dim == 1
               M = sparse((-3*q(1)./2.*eye(Nx)+D));
       else
               M = sparse((-3*q(1)./2.*eye(Nx*Nx)+D));
       end
    end

    for iteration = progress(1:T/dt)
% % % % % % % % % %         BDF1
        if order == 1
            
            if iteration == 1
                
                b = source-q(1).*history_arr(:,end);
                
                root = M\b;
                
            else
                
                diff_b = diff(history_arr, 1, 2);

                q_r = q(iteration:-1:2);

                b = source+sum(diff_b.*q_r, 2)-q(1).*history_arr(:,end);
%                 for time = 1:iteration-1
%                      
%                      b = b+q(time+1).*(history_arr(:,end-(time-1))-history_arr(:,end-(time-1)-1));     
%                      
%                 end
                root  = M\b;
            end
% % % % % % % % % %         BDF2
        elseif order == 2
            

            if iteration == 1
                
                b = source+q(1).*(history_arr(:,end-1)-4.*history_arr(:,end))./2;
                root  = M\b;
                
            else
%                 diff_b = (3.*history_arr(:,3:end)-4.*history_arr(:,2:end-1)+history_arr(:,1:end-2))./2;
                diff_b_tmp = diff(history_arr, 1, 2);
                
                diff_b = (3.*diff_b_tmp(:, 2:end)-diff_b_tmp(:, 1:end-1))./2;
                q_r = q(iteration:-1:2);

                b = source+sum(diff_b.*q_r, 2)-q(1).*(4.*history_arr(:,end)-history_arr(:,end-1))./2;
                root = M\(b);
            end
        end

        if sum(abs(root)>1)
            disp('rrr')
            break
        else
            history_arr = [history_arr, root(:)];
        end
    end
%     history_arr = history_arr(:,end);
end 