function K_est = Sindy_rossler(x,x_dot,r)
    
    lambda = 0.002; 
    
    x_dot = x_dot(:,1:3:end);
    N = size(x_dot,2);
    M = size(x_dot,1);
    K_est = zeros(N,N);
    
    for i = 1:N
        A = ones(M,1);
        
        % Add the basis function expension for the self dyanmics 
        for j = 1:r-1
            A = [A x(:,3*(i-1)+1).^j x(:,3*(i-1)+2).^j x(:,3*i).^j];
        end
        
      
        x_coupling = x(:,1:3:end);
        
        % Add coupling dynamics
        for j = 1:r   
            for k = 0:r-1
                if(i==1)
                    A = [A (x_coupling(:,i).^k).*(x_coupling(:,2:end).^j)];
                else
                    A = [A (x_coupling(:,i).^k).*(x_coupling(:,1:i-1).^j)];
                    A = [A (x_coupling(:,i).^k).*(x_coupling(:,i+1:end).^j)];
                end
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
                K_est(i,j) = norm(Coeff(3*(r-1)+j+1:N-1:end),'fro');
            else
                K_est(i,j) = norm(Coeff(3*(r-1)+j:N-1:end),'fro');
            end
        end
    end
    
end