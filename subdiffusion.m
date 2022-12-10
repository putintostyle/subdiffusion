function history_arr = subdiffusion(order, Nx, D, init, alpha, T, dt, dim)
    
    q = q_weight(T/dt, dt, alpha);
    if dim == 1
        history_arr = [zeros(Nx,1),zeros(Nx,1)];
    else
        history_arr = [zeros(Nx*Nx,1),zeros(Nx*Nx,1)];
    end

    source = -D*init;
    for iteration = progress(1:T/dt)
% % % % % % % % % %         BDF1
        if order == 1
            b = -1*q(1).*history_arr(:,end);
            if iteration == 1
                if dim == 1
                    root  = linsolve((-q(1).*eye(Nx)+D),(b-source));
                else
                    root  = linsolve((-q(1).*eye(Nx*Nx)+D),(b-source));
                end
            else
                for time = 1:iteration-1
                     
                     b = b+q(time+1).*(history_arr(:,end-(time-1))-history_arr(:,end-(time-1)-1));     
                     if dim == 1
                        root  = linsolve((-q(1).*eye(Nx)+D),(b-source));
                     else
                        root  = linsolve((-q(1).*eye(Nx*Nx)+D),(b-source));
                     end
                end
            end
% % % % % % % % % %         BDF2
        elseif order == 2
            b = q(1).*(history_arr(:,end-1)-4.*history_arr(:,end))./2;
            if iteration == 1
                if dim == 1
                    root  = linsolve((-3*q(1)/2.*eye(Nx)+D),(b-source));
                else
                    root  = linsolve((-3*q(1)/2.*eye(Nx*Nx)+D),(b-source));
                end
            else
                for time = 1:iteration-1
                     
                     b = b+q(time+1)*(3.*history_arr(:,end-(time-1))-4.*history_arr(:,end-(time-1)-1)+history_arr(:,end-(time-1)-2))/2;     
                     if dim == 1
                        root  = linsolve((-3*q(1)./2.*eye(Nx)+D),(b-source));
                     else
                        root  = linsolve((-3*q(1)./2.*eye(Nx*Nx)+D),(b-source));
                     end
                end
            end
        end

        if sum(abs(root)>1)
            root
            break
        else
            history_arr = [history_arr, root(:)];
        end
    end
    history_arr = history_arr(:,end);
end 