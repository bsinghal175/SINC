function K_est = LASSO_kuramoto(x,x_dot,r)
    

    %Infer connectivity using LASSO regression

    N = size(x,2);
    M = size(x,1);
    K_est = zeros(N,N);
    
    for i = 1:N
        A = ones(M,2*N*r+1);
        x_diff = x(:,:) - repelem(x(:,i),1,N);
        data_A = kron(linspace(1,r,r),x_diff);
        A(:,2:2*N*(r)+1) = [sin(data_A) cos(data_A)];
        A(:,1+i:N:end) = [];
        
        A = A(:,2:end);
        [z,FitInfo] = lasso(A,x_dot(:,i),'CV',5,...
            'Options',statset('UseParallel',true));
        idxLambda1SE = FitInfo.Index1SE;
        z = z(:,idxLambda1SE); 
        z = [0;z];
          
        for j = 1:N
            if i == j
                continue
            elseif j < i
                coupling_coeff = z(j+1:N-1:end);
                K_est(i,j) = norm(coupling_coeff,'fro');
            else
                coupling_coeff = z(j:N-1:end);
                K_est(i,j) = norm(coupling_coeff,'fro');
            end
        end  
    end


    
end