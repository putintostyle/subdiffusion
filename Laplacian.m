function D = Laplacian(nx, ny, h_step ,dim, method)
    if method == 0 % fdm
        if dim == 1
            I = eye(nx);
            row = 2:nx;
            col = 1:nx-1;
            data = ones(1,nx-1);
            off_diag = full(sparse(row, col, data, nx, nx));
            off_diag = off_diag+off_diag';
            D = off_diag-2.*(eye(nx));
            D = D./h_step./h_step;
        elseif dim == 2
            I = eye(nx);
            row = 2:nx;
            col = 1:nx-1;
            data = ones(1,nx-1);
            
            off_diag = full(sparse(row, col, data, nx, nx));
            off_diag = off_diag+off_diag';
            S = off_diag-4.*(eye(nx));
            D = kron(I, S)+kron(off_diag, I);
            D = D./h_step./h_step;
    
        end

    elseif method == 1 % fdm
        if dim == 1
            I = eye(nx);
            row = 2:nx;
            col = 1:nx-1;
            data = ones(1,nx-1);
            off_diag = full(sparse(row, col, data, nx, nx));
            off_diag = off_diag+off_diag';
            D = off_diag-2.*(eye(nx));
            D = D./h_step./h_step;
        elseif dim== 2
            row = 2:nx^2;
            col = 1:nx^2-1;
            row_2 = 4:nx^2;
            col_2 = 1:nx^2-3;
            data = ones(1,nx^2-1);
            data_2 = ones(1,nx^2-3); 
            off_diag = full(sparse(row, col, data, nx^2, nx^2));
            off_diag = off_diag+off_diag';
            off_diag_2 = full(sparse(row_2, col_2, data_2, nx^2, nx^2));
            off_diag_2 = off_diag_2+off_diag_2';
            D = off_diag_2+off_diag-4.*eye(nx^2);
    
            D = D./h_step./h_step;
        end

            
    end


         

end