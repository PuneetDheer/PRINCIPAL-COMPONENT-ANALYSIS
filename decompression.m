%Function to Get back the original data or almost original data
%%
function [ Orig_data ] = decompression( choice, Transformed_data, Eigen_vector_col_mat, Original_mean, Original_std )


    if (choice == 'demean')
            Orig_data = bsxfun( @plus,( Transformed_data * Eigen_vector_col_mat' ), Original_mean ); 

    elseif (choice =='zscore')
            Orig_data = bsxfun( @times,( Transformed_data * Eigen_vector_col_mat' ), Original_std ); 
            Orig_data = bsxfun( @plus, Orig_data, Original_mean );

    end

end

