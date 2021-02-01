function [Ind_depth, F_m_col, Contact_point, Max_displacement, Contact_point_row, diff] = Contact_point( B_Extend,number_curves )
%Calcuation of Contact point on extend curve
% Contact point is determinded by fitting a differential to every data point along the extend curve
% When the differential of the next point is 1.5x > previous point, this is the CP 

    for i=1:number_curves
        F_m(i) = max(B_Extend{i}(:,2)); 
        F_m_col(i) = find(B_Extend{i}(:,2)==F_m(i));
        Max_displacement(i) = B_Extend{i}(F_m_col(i),1);
        for j = 2:length(B_Extend{i})
            diff{i}(j) = (B_Extend{i}(j,2) - B_Extend{i}(j-1,2)) / ((B_Extend{i}(j-1,1)) - (B_Extend{i}(j,1)));
            if (diff{i}(j)) > (1.5*diff{i}(j-1));  
                Contact_point(i) = B_Extend{i}(j-1,1);
                Contact_point_row(i) = j-1;
            end
        end
        
        Ind_depth(i) =  Contact_point(i) - Max_displacement(i);
    end
end

