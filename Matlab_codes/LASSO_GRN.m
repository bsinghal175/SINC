function K_est = LASSO_GRN(x,x_dot,r)
    
    N = size(x,2);
    M = size(x,1);
    K_est = zeros(N,N);
    
    for i = 1:N
        A = ones(M,1);
        
        % Add the basis function expension for the self dyanmics 
        for j = 1:r
            A = [A x(:,i).^j];
        end
        
        % Add coupling dynamics
        for j = 1:r
            if(i==1)
                A = [A x(:,2:end).^j];
            else
                A = [A x(:,1:i-1).^j];
                A = [A x(:,i+1:end).^j];
            end
        end
        
        A = A(:,2:end);
        [z,FitInfo] = lasso(A,x_dot(:,i),'CV',5,...
            'Options',statset('UseParallel',true));
        idxLambda1SE = FitInfo.Index1SE;
        z = z(:,idxLambda1SE); 
        Coeff = [0;z];
        
        
        %Estimate the coupling matrix
        for j = 1:N
            if(i==j)
                continue
             elseif j < i
                K_est(i,j) = norm(Coeff(r+j+1:N-1:end),'fro');
            else
                K_est(i,j) = norm(Coeff(r+j:N-1:end),'fro');
            end
        end
    end
    
end