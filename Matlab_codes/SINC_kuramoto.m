function results = SINC_kuramoto(x,x_dot,r,epsilon,Gamma)
    
    % x : time-series data of size MxN. 
    % x_dot : Derivative of the data
    % r : Number of Fourier terms used for approximation
    % epsilon: Convergence criteria 
    % Gamma : number of network iterations

    N = size(x,2);
    K_est = zeros(N,N);
    
    %Initialize and normalize the Fourier Coefficients
    Lambda_initial = ones(2*r,1);
    Lambda_initial = (1/norm(Lambda_initial) )*Lambda_initial;
    
    Iteration_number = zeros(N,Gamma);
    Lambda_hat = Lambda_initial;
    Lambda_each_node = zeros(N*Gamma,2*r);
    %B_each_node = zeros(N*Gamma,N-1);
    %A_each_node = zeros(N*Gamma,r);
    for network_iteration = 1:Gamma
        % i denotes the node number  
        for i = 1:N
          
            F = Estimate_F_matrix(x(:,i),r);
            Z = Estimate_Z_matrix(x,r,i);
     
            
            
            % Estimate the connection strength(B) and the self-dynamics(A)
            for iter = 1:100
                [B_hat,A_hat] = Estimate_topology(x_dot(:,i),F,Z,Lambda_hat,N);
                Lambda_hat_updated = Estimate_lambda(x_dot(:,i),F,Z,B_hat,A_hat,N,r);
                
                k = norm(Lambda_hat_updated);
                Lambda_hat_updated = Lambda_hat_updated/k;
                B_hat = k*B_hat;
                
                if(norm(Lambda_hat_updated-Lambda_hat)<epsilon || iter == 100)
                    Lambda_hat = Lambda_hat_updated;
                    Iteration_number(i,network_iteration) = iter;
      
                    
                    % Store the values
                    Lambda_each_node((network_iteration-1)*N+i,:) = Lambda_hat';
                    break;
                end
                Lambda_hat = Lambda_hat_updated;
            end
        end
    end
    
    Lambda = mean(Lambda_each_node((Gamma-1)*N+1:Gamma*N,:));
    for i = 1:N
        F = Estimate_F_matrix(x(:,i),r);
        Z = Estimate_Z_matrix(x,r,i);
        [B_hat(:,i),A_hat(:,i)] = Estimate_topology(x_dot(:,i),F,Z,Lambda',N);
    end
             
    % Estimate the coupling matrix
    for i = 1:N
        for j= 1:N
            if(i<j)
                K_est(i,j) = B_hat(j-1,i);
            elseif(i>j)
                K_est(i,j) = B_hat(j,i);
            end
        end
    end
           
    results.K_est = K_est;
    results.A_est = A_hat;
    results.Lambda = Lambda;
    results.iteration = Iteration_number;
    results.Lambda_iteration = Lambda_each_node;


end


%% Functions used in the main file
function Z = Estimate_Z_matrix(x,r,node_index)
    M = size(x,1);
    N = size(x,2);
    
    Z = zeros(M*(N-1),2*r);
    
    x_temp = x;
    x_temp(:,node_index) = [];
    for t = 1:M
        Z_i = zeros(N-1,2*r);
        
            for k = 1:r
                Diff = x(t,node_index)-x_temp(t,:);
                Z_i(:,2*k-1:2*k) = [sin(k*Diff)' cos(k*Diff)'];
            end
        
        Z((t-1)*(N-1)+1:t*(N-1),:) = Z_i;
    end
    
end

function F = Estimate_F_matrix(x,r)
    
    M = size(x,1);
    F = ones(M,1);
%     for j = 1:r
%         F(:,j) = x.^j;
%     end 
end


function [B_hat,A_hat] = Estimate_topology(x_dot,F,Z,Lambda,N)
    M = size(x_dot,1);
    X = zeros(M,N-1);
   
    for i = 1:M
        X(i,:) = (Z((i-1)*(N-1)+1:i*(N-1),:)*Lambda)';
    end
    
    Temp = eye(M)-(F*pinv(F));
    
    B_hat = lsqnonneg(Temp*X,Temp*x_dot);
    A_hat = pinv(F)*(x_dot-X*B_hat);
end

function Lambda_hat = Estimate_lambda(x_dot,F,Z,B_hat,A_hat,N,r)
    M = size(x_dot,1);
    X = zeros(M,2*r);
    
    for i = 1:M
        X(i,:) = B_hat'*Z((i-1)*(N-1)+1:i*(N-1),:);
    end
    
    Lambda_hat = pinv(X)*(x_dot-F*A_hat);
end