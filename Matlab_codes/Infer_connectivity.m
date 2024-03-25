%% This code inferes the network connectivity structure using the SINC 
%  approach. 

% Inputs: 
%       1. Model_type - 1 for Kuramoto oscillators; 2 for gene-regulatory
%       networks & 3 for Rossler oscillators 

% The functions LASSO_"model name" and Sindy_"model name" infer
% connectivity using LASSO or Sindy. 


Model_type = 3;  %change this to infer other network models 
switch(Model_type)
    
     
    %---------------------------------------------------------------------
    case 1  %Infer the connectivity of Kuramoto oscillators
        
        %Algorithm hyperparameters
        Gamma = 2;
        Epsilon = 10^-3;
        r = 2;
          
        X = readmatrix('..\Data\Kuramoto\Time_series.csv');
        X_dot = readmatrix('..\Data\Kuramoto\Time_series_derivative.csv');
        Adjacency_matrix =  readmatrix('..\Data\Kuramoto\Adjacency_matrix.csv');
        Adjacency_matrix_binary = Adjacency_matrix;
        Adjacency_matrix_binary(Adjacency_matrix_binary>0) = 1;

        %Infer networks for 6 different data lengths 
        for i = 1:6
            results(i) = SINC_kuramoto(X(1:50*i,:),...
                X_dot(1:50*i,:),r,Epsilon,Gamma);
            [~,~,~,AUC(i)] = perfcurve...
                (Adjacency_matrix_binary(:),results(i).K_est(:),1);   
            
            % To infer connectivity using LASSO 
            K_est = LASSO_kuramoto(X(1:50*i,:),...
                X_dot(1:50*i,:),r);
            [~,~,~,AUC_lasso(i)] = perfcurve...
                (Adjacency_matrix_binary(:),K_est(:),1);  
            
            % To infer connectivity using Sindy 
            K_est = Sindy_kuramoto(X(1:50*i,:),...
                X_dot(1:50*i,:),r);
            [~,~,~,AUC_Sindy(i)] = perfcurve...
                (Adjacency_matrix_binary(:),K_est(:),1);
            
        end
      
    
        
    %---------------------------------------------------------------------
    case 2  %Infer the connectivity of Gene Regulatory networks  
        %Algorithm hyperparameters
        Gamma = 2;
        Epsilon = 10^-3;
        r = 3;
        
        X = readmatrix('..\Data\Gene_regulatory\Time_series.csv');
        X_dot = readmatrix('..\Data\Gene_regulatory\Time_series_derivative.csv');
        Adjacency_matrix =  readmatrix('..\Data\Gene_regulatory\Adjacency_matrix.csv');
        Adjacency_matrix_binary = Adjacency_matrix;
        Adjacency_matrix_binary(Adjacency_matrix_binary>0) = 1;
        
        %Infer networks for 8 different data lengths 
        count = 1;
        for i = 1:2:15
            results(count) = SINC_GRN(X(1:40*i,:),...
                X_dot(1:40*i,:),r,Epsilon,Gamma);
            [~,~,~,AUC(count)] = perfcurve...
                (Adjacency_matrix_binary(:),results(count).K_est(:),1);  
            
            
            
            % To infer connectivity using LASSO 
            K_est = LASSO_GRN(X(1:40*i,:),...
                X_dot(1:40*i,:),r);
            [~,~,~,AUC_lasso(count)] = perfcurve...
                (Adjacency_matrix_binary(:),K_est(:),1);  
            
            % To infer connectivity using Sindy 
            K_est = Sindy_GRN(X(1:40*i,:),...
                X_dot(1:40*i,:),r);
            [~,~,~,AUC_Sindy(count)] = perfcurve...
                (Adjacency_matrix_binary(:),K_est(:),1);
            
            count = count+1;
        end
     
        
    %---------------------------------------------------------------------
    case 3  %Infer the connectivity of Rossler oscillator networks
        %Algorithm hyperparameters
        Gamma = 2;
        Epsilon = 10^-3;
        r = 3;
        
        X = readmatrix('..\Data\Rossler\Time_series.csv');
        X_dot = readmatrix('..\Data\Rossler\Time_series_derivative.csv');
        Adjacency_matrix =  readmatrix('..\Data\Rossler\Adjacency_matrix.csv');
        Adjacency_matrix_binary = Adjacency_matrix;
        Adjacency_matrix_binary(Adjacency_matrix_binary>0) = 1;
        
        %Infer networks for 6 different data lengths 
        for i = 1:6
            results(i) = SINC_rossler(X(1:50*i,:),...
                X_dot(1:50*i,:),r,Epsilon,Gamma);
            [~,~,~,AUC(i)] = perfcurve...
                (Adjacency_matrix_binary(:),results(i).K_est(:),1);   
            
            
            % To infer connectivity using LASSO 
            K_est = LASSO_rossler(X(1:50*i,:),...
                X_dot(1:50*i,:),r);
            [~,~,~,AUC_lasso(i)] = perfcurve...
                (Adjacency_matrix_binary(:),K_est(:),1);  
            
            % To infer connectivity using Sindy 
            K_est = Sindy_rossler(X(1:50*i,:),...
                X_dot(1:50*i,:),r);
            [~,~,~,AUC_Sindy(i)] = perfcurve...
                (Adjacency_matrix_binary(:),K_est(:),1);
        end
        
end