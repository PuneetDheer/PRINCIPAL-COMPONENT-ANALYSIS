% PRINCIPAL COMPONENT ANALYSIS (PCA)
% Coded By: Puneet Dheer
% Date:15-NOV-2018
%
% INPUT:
% choice = for standardizing the data ['demean' or 'zscore']
% Column_data = Column wise data (column represents dimension)
%
% OUTPUT:
% Reduced_Transformed_data = compressed data in reduced dimension
% All_Transformed_data = compressed data in all dimension
% Reconst_Orig_data = reconstructed data from All_Transformed_data
% Cov_mat = Covariance matrix of standardized input data
% Eigen_vector_col_mat = eigen vector of  Covariance matrix in sorted order (PC1--to--PCn) w.r.t. maximum VARiance or Eigen_value
% Eigen_value = eigen vector of  Covariance matrix in descending order 
% VARiance = latent or variance of All_Transformed_data  
% PC_Confidence = Percentage of Variation
% {Note: Eigen_value will always be equal to VARiance, so both the parameter can be used to identity the highest variance. Here, VARiance was used.}
%%

function [ Reduced_Transformed_data, All_Transformed_data, Reconst_Orig_data, Cov_mat, Eigen_vector_col_mat, Eigen_value, VARiance, PC_Confidence ] = PCAnalysis( choice, Column_Data )

%{1} Original data consists of 'm' data points(in rows) and 'n' Dimensional data(in columns) 
%{2} PCA Goal is to reduce the 'n' Dimension into 'k' (k<=n) while preserving most of the information about original data.


    %Step 1: DEmean  or zscore the data (always) --> centering the data to the origin
    if (choice == 'demean')
        X = bsxfun(@minus,Column_Data,mean(Column_Data));
    elseif (choice == 'zscore')
        X=zscore(Column_Data);
    end

    Original_mean = mean(Column_Data);
    Original_std = std(Column_Data);
    Original_dimension = size(Column_Data,2); 

    %Step 2: Calculate the Covariance Matrix
    Cov_mat = cov(X);

    %Step 3: Calculate Eigen value and corresponding eigen vector
    [Eigen_vector_col_mat, Diag_Eigen_value] = eig(Cov_mat); % k eigen vector (k == n)
    Eigen_value = diag(Diag_Eigen_value)';
%     [~,ind] = max(abs(Eigen_vector_col_mat)); 
%     set_sign = sign(Eigen_vector_col_mat(sub2ind(size(Eigen_vector_col_mat),[ind],[1:Original_dimension])));
%     Eigen_vector_col_mat = bsxfun(@times, Eigen_vector_col_mat, set_sign); %mutiply with constant -1 only flip the vector

    %Step 4: Multiply the X data with eigen vector (transforming the data into 'k' independent dimension)
    All_Transformed_data = X * Eigen_vector_col_mat; %Compressing the data X into 'k' independent eigen vector dimension

    %Step 5: Calculate the Variance (latent) of all_Transformed_data (it will be equal to the eigen value of Covariance Matrix)
    VARiance = var(All_Transformed_data);

    %Step 6: Calculate confidence in percentage (it will help in selecting the most prominent dimension, the higher will be better)
    Confidence = (VARiance./sum(VARiance))*100;
    [PC_Confidence, PC_eigenvector_order] = sort(Confidence,'descend');


    %Final results will be in sorted order from PC1...to...PCn, based on the sorted value of  VARiance or Eigen_value (from highest to lowest eigen value)
    All_Transformed_data = All_Transformed_data(:,PC_eigenvector_order);
    Eigen_vector_col_mat = Eigen_vector_col_mat(:,PC_eigenvector_order);
    Eigen_value = Eigen_value(:,PC_eigenvector_order);
    VARiance = VARiance(:,PC_eigenvector_order);
    
    figure
    h1 = bar(PC_Confidence);
    xlabel('Principal Components (Dimensions)')
    ylabel('Percentage of Variance Explained')
    hold on
    h2 = plot(cumsum(PC_Confidence),'r-o','MarkerSize',8,'LineWidth',2,'MarkerEdgeColor','g','MarkerFaceColor','r');
    legend([h1 h2],{'Individual Explained Variance', 'Cumulative Explained Variance'})
    title('Scree Plot')
    set(gca,'FontSize',11, 'FontWeight','bold')
    grid on
    axis tight

    msg = sprintf('Enter the no. of reduced dimension in range : = %d to %d \n', 1, Original_dimension );
    k = input(msg);
    Reduced_Transformed_data = All_Transformed_data(:,1:k);

    %Step 7: [Reconstruction] => Get back the {original data (if (k == n))}  or {almost the original data (if (k < n))}  (it will loose some information) after deciding 
    %which eigen vector will be most suitable without loosing so much information.
    Reconst_Orig_data = decompression( choice, All_Transformed_data(:,1:k), Eigen_vector_col_mat(:,1:k), Original_mean, Original_std );


end
