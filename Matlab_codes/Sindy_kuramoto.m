function K_est = Sindy_kuramoto(x,x_dot,r)
    
    %Infer connectivity using Sindy algorithm
    
    lambda = 0.002; 
    N = size(x,2);
    M = size(x,1);
    K_est = zeros(N,N);
 
    for i = 1:N
        A = ones(M,2*N*r+1);
        x_diff = x(:,:) - repelem(x(:,i),1,N);
        data_A = kron(linspace(1,r,r),x_diff);
        A(:,2:2*N*(r)+1) = [sin(data_A) cos(data_A)];
        A(:,1+i:N:end) = [];
        
        dXdt = x_dot(:,i);
        z = A\dXdt;  % initial guess: Least-squares

        % lambda is our sparsification knob.
        for k=1:20
            smallinds = (abs(z)<lambda);   % find small coefficients
            z(smallinds)=0;                % and threshold
            biginds = ~smallinds;
            % Regress dynamics onto remaining terms to find sparse Xi
            z(biginds) = A(:,biginds)\dXdt; 
            
        end
        
        
  
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