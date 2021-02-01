%% Clear previous workspace 

clc 
clear all 

%% Setting file parameters  

number_curves=10; % Choose number of curves to import
percentage=0.25; % Chose percentage of full retract curve to take to calculate slope (25% is a good place to start for this)
startingrow = 79; % Starting row in the txt file to import the data 
Folder = 'Raw_Curves'; % Folder where the force curve text files are saved 


%Columns to Import

vertical_tip_position =1; % This should be vertical tip position (m) (X-axis)
vertical_deflection = 2; % This is vertical deflection  (N)(Y-axis)

% Parameters to set for Oliver & Pharr model 

R= 5e-6; % Bead radius (meters)
v=0.5; % Poisson's ratio 
    
 %% Load extend data 
 
 for i =1:number_curves; 
    ivalue=int2str(i);
    f_Extend=fullfile(Folder,['extend', ivalue, '.txt']); % Load extend curves
    A_Extend{i}= importdata(f_Extend,' ',startingrow); %Import data from text file
    B_Extend{i}=[A_Extend{1,i}.data(:,vertical_tip_position),A_Extend{1,i}.data(:,vertical_deflection)]; % Save Vertical Tip Position and Vertical Deflection into a new matrix
      
    
 hold on
 
 figure(1)
 plot(B_Extend{i}(:,1),B_Extend{i}(:,2)); %% Plot the extend curves
 title('Extend curves');
 xlabel('Vertical tip position (m)') 
 ylabel('Vertical deflection (N)')
 hold on 
 
%Find Max and Min Force 
   
maxforce_extend(i)=max(B_Extend{i}(:,2)); %% Find Pmax (Max Force) 
maxforce_extend_row(i)=find(B_Extend{i}(:,2)==max(B_Extend{i}(:,2))); %Find Row Pmax is located
   
minforce_extend(i)= min(B_Extend{i}(:,2)); %% Find min force (Fmin) 
minforce_extend_row(i)=find(B_Extend{i}(:,2)==min(B_Extend{i}(:,2))); %Find Row Fmin is located
    
 end
 
% Run Contact point function 

[Ind_depth, F_m_col, Contact_point, Max_displacement, Contact_point_row, diff] = Contact_point(B_Extend,number_curves);
 

    
%% Load retract data


for i =1:number_curves; 


    ivalue=int2str(i);
    f_Retract=fullfile(Folder,['retract', ivalue, '.txt']); % Load retract curves 
    A_Retract{i}= importdata(f_Retract,' ',startingrow); %Import data from text file
    B_Retract{i}=[A_Retract{1,i}.data(:,vertical_tip_position),A_Retract{1,i}.data(:,vertical_deflection)]; % Save Vertical Tip Position and Vertical Deflection to a new matrix
        
 hold on
 
 figure(2)
 plot((B_Retract{i}(:,1)),(B_Retract{i}(:,2))); %% Plotting the retract curve 
 title('Retract curves');
 xlabel('Vertical tip position (m)') 
 ylabel('Vertical deflection (N)')
 hold on
   

 %Find Max and Min Force 

maxforce_retract(i)=max(B_Retract{i}(:,2)); %% Find Fmax 
maxforce_retract_row(i)=find(B_Retract{i}(:,2)==max(B_Retract{i}(:,2))); %Find row Fmax (Max Force) is located 

minforce_retract(i)= min(B_Retract{i}(:,2)); %% Find fmin 
minforce_retract_row(i)=find(B_Retract{i}(:,2)==min(B_Retract{i}(:,2))); %Find Row Fmin is located
 

TipPosition_Fmin(i)=B_Retract{i}((minforce_retract_row(i)),1); % Find the vertical tip position (um) for Fmin vertical deflection
 
TipPosition_Fmax(i)=min(B_Retract{i}(maxforce_retract_row(i),(1))); %Find the vertical tip position (um) for Fmax vertical defelction 


% Finding the Retract indentation dept: Fmax - Fmin
        
Ind_depth_retract(i)=(TipPosition_Fmin(i) - TipPosition_Fmax(i)); % Indentation range: Fmax to Fmin
 
indent_range_rows(i)=(minforce_retract_row(i) - maxforce_retract_row(i)); %% Indentation range based on row number 

%Comment: Indent range looks at Fmax - Fmin - so a small amount of adhesion is neccesary to be negative so that you can find the slope
      
Pmax(i)=maxforce_extend(i); % Use Pmax from extend curve


 %% Calculating the Young's modulus 
 
 %Select upper percentage of curve to fit (Sectioned Curve)
   
   gradbegin(i)=round(percentage*indent_range_rows(i)); % Chosen percentage of where to start (row number) for gradient
    
   
   Force_SectionedCurve(i)=(B_Retract{i}(gradbegin(i),2)); %% Vertical deflection at gradient begin row 
    
   SectionedCurve_TipEnd(i)=(B_Retract{i}(gradbegin(i),1)); % Vertical tip position at gradient begin row
   
   SectionedCurve_TipStart(i)=TipPosition_Fmax(i);  % Tip start at Fmax in retract curve
  

 %Fitting the Curve - Polynomial fit (1st order)
       
poly_fit{i}=[polyfit((B_Retract{i}(1:gradbegin(i),1)),(B_Retract{i}(1:gradbegin(i),2)),1)];  % This does a poly fit for the sectioned curve


grad_fitted(i)=abs(poly_fit{i}(1,1)); % Take the gradient (m) value from polyfit and make it possitive
    

%Manually calculate the gradient for the sectioned curve - Check point 
 
SectionedCurve_gradient(i)=((maxforce_retract(i)-Force_SectionedCurve(i))/(SectionedCurve_TipEnd(i)-SectionedCurve_TipStart(i)));  %This divides the absolute of the vertical deflection range / indentation depth range (e.g. the gradient)
    

%Manually calculate the gradient for the full retract cuve (Fmax - Fmin) 
 
FullCurve_gradient(i)=abs(Pmax(i)/Ind_depth_retract(i)); 
   

%Plot retract curve (Fmin - Fmax)

figure(3) 
plot((B_Retract{i}((maxforce_retract_row(1):minforce_retract_row(i)),1)),(B_Retract{i}((maxforce_retract_row(i):minforce_retract_row(i)),2)));
title('Full retract curve (Fmax - Fmin)');
xlabel('Vertical tip position (m)') 
 ylabel('Vertical deflection (N)')
hold on;
    
%Plot sectioned retract curve 

 figure(4)
 plot((B_Retract{i}((1:gradbegin(i)),1)),(B_Retract{i}((1:gradbegin(i)),2))); %This plots the slope for the sectioned gradient (Percentage) 
 title('Retract sectioned');
 xlabel('Vertical tip position (m)') 
 ylabel('Vertical deflection (N)')
 hold on;
    
   
%% Calculating parameters for the Oliver & Pharr model: 

    
grad_used(i)=grad_fitted(i); % This takes the gradient from the polyfit 
    
epsilon=0.75; % This depends on the indenter geometry - 0.75 for a paraboloid of revolution  
 
    
h_max_retract(i)=abs(Ind_depth_retract(i)); %% This is length (um)from vertical tip position Fmax - Fmin  


h_max (i) = abs(Ind_depth(i)); % Hmax (maximum displacement) - Calculated from the extend curve
 

h_s(i)=epsilon*(Pmax(i)/grad_used(i)); %This calculates h-s from: Epsilon * (Pmax (Max force) / S (gradient fitted))
    

h_c(i)=abs(h_max(i)-h_s(i)); %%  Calculate h-c from h_max and h_s
    
        
A_hc(i)=(pi*((2*R*h_c(i)) - (h_c(i)^2))); % Calculates contact area for tip


num_eq(i)=grad_used(i)*(sqrt(pi)); % This is the numerator of the Oliver and Pharr equation (Supplementary note 1)

den_eq(i)=2*sqrt(A_hc(i)); %This is the denominator of the Oliver and Pharr equation (Supplementary note 1)
  
Eop_row(i)=num_eq(i)/den_eq(i); % Reduced modulus  -Does not take into account Poisson's ratio
    
  
    Eop_column=Eop_row';
    EopkPa_row(i)=Eop_row(i)/1000; % Young's modulus in kPa
    EopkPa_column=EopkPa_row'; % Young's modulus in kPa
    

 %Calculate Young's modulus - based on Poisson's ratio
     
  Youngs_Modulus(i) = Eop_column(i)*(1-v^2); 
 
  % This assumes no deformation of the indentor ! 
  

end 
 

%% Hisogram of Indentation depth and Full retract length 
figure(5)
histogram(Ind_depth, 'Edgecolor', 'b');
hold on;
histogram(Ind_depth_retract, 'Edgecolor', 'r');
grid on;
% Put up legend.
legend1 = sprintf('Indentation depth', mean(Ind_depth));
legend2 = sprintf('Full retract length', mean(Ind_depth_retract));
legend({legend1, legend2,});
ylabel('Frequency') 
xlabel('Indentation depth/length (m)')


%% Save data to excel files 

File = 'Youngs modulus.xls'
xlswrite (File,Youngs_Modulus') % Young's modulus in pascal

File = 'Reduced modulus.xls'
xlswrite (File,Eop_column) % Reduced modulus in pascal





