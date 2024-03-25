function K_est = Sindy_GRN(x,x_dot,r)
    
    lambda = 0.002; 
    
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
        
        Coeff = z;
        
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