% Solves temperatures within model.
% Script File: temperatureSolverInitialisation.m
% Author: Stephen Blowers  -  S.Blowers@ed.ac.uk
% Date Modified: 08/02/2017
% Description: Solves temperatures using initialised options and
% parameters. Runs the script temperatureSolverSetUp or 
% temperatureSolverSetUp_Pennes depending on what options are chosen to
% create a sparse matrix that is then solved as a linear set of equations.



disp('%%%%%%% Solving Temperatures %%%%%%%%%%%%%%%%%%%%%%') % Display.

%%%%%% Using Graetz Number %%%%%%%%%%%%%%%%%%%%%%

if Option_GraetzNumber && ~Option_PennesOnly
    Gz1 = Rho_b*Cp_b*(max(abs(FdotArt),[],2)/Rho_b./Aavg1).*Davg1.^2/Kc_b./L1;
    Gz2 = Rho_b*Cp_b*(min(abs(FdotVein),[],2)/Rho_b./Aavg2).*Davg2.^2/Kc_b./L2;
    Nu1 = Nu+0.155*exp(1.58*log10(Gz1));
    Nu2 = Nu+0.155*exp(1.58*log10(Gz2));
    
    Nu1(Nu1>20) = 20;
    Nu2(Nu2>20) = 20; % Capped as the equation is only valid for Gz<1000
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Vessel Inter-Domain Heat Transfer Terms %%
if ~Option_PennesOnly
    Beta1 = zeros(size(Vessel1,1),1);
    for n = 2:size(Vessel1,1)
        if size(Vessel2Volume1{n},1)~=0
            if Option_GraetzNumber
                Beta1(n) = (Nu1(n)*Kc_b./Davg1(n)).*AL1(n) / size(Vessel2Volume1{n},1);
            else
                Beta1(n) = (Nu*Kc_b./Davg1(n)).*AL1(n) / size(Vessel2Volume1{n},1);
            end
        else
            Beta1(n) = 0;
        end
    end
    % Setting up inter-domain heat transfer terms for each arterial line
    % segment.
    
    Beta2 = zeros(size(Vessel2,1),1);
    for n = 2:size(Vessel2,1)
        if size(Vessel2Volume2{n},1)~=0
            if Option_GraetzNumber
                Beta2(n) = (Nu2(n)*Kc_b./Davg2(n)).*AL2(n) / size(Vessel2Volume2{n},1);
            else
                Beta2(n) = (Nu*Kc_b./Davg2(n)).*AL2(n) / size(Vessel2Volume2{n},1);
            end
        else
            Beta2(n)=0;
        end
    end
    % Setting up inter-domain heat transfer terms for each venous line
    % segment.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Fixing LimitedCoolingHeight if Needed %%%%
if Option_LimitedCooling
    if strcmp(Option_DomainCoarsening, 'two') || strcmp(Option_DomainCoarsening, 'six')
        LimitedCoolingHeight = round(LimitedCoolingHeight/2);
    end
    if strcmp(Option_DomainCoarsening, 'three') || strcmp(Option_DomainCoarsening, 'six')
        LimitedCoolingHeight = round(LimitedCoolingHeight/3);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Interactions %%%%%%%%%%%%%%%%%%%
if Option_PennesOnly
    temperatureSolverSetUp_PennesOnly
else
    temperatureSolverSetUp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
%%%%%% Solving Linear System %%%%%%%%%%%%%%%%%%%%
disp('Starting Linear Solver') % Display.
disp(['Solving Matrix Inversion size: ' num2str(size(T_Solve,1)) ' by ' num2str(size(T_Solve,2))]) % Display.
tic % Start timing.

T_N = T_Solve\D; % Solve the linear system.

toc % Finish timing.
disp('Matrix Inversion Completed') % Display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Reshaping Resutls from Solver %%%%%%%%%%%%
if Option_PennesOnly
    Tt = zeros(size(DomTot));
    Tt(RowConvert) = T_N(1:NumDomTot);
    
    Tt(~DomTot) = NaN;
    % Setting all values outside of the domain to NaN. This facilitates
    % visualisation of data with 'planecut'.
else
    Tt = zeros(size(DomTot));
    Tt(RowConvert) = T_N(1:NumDomTot);
    
    if Beta23 == Inf
        Tb = Tt;
    else
        Tb = T_N(NumDomTot+1:NumDomRows);
    end
    
    T_Art = T_N(NumDomRows+1:NumDomRows+VesselRow);
    T_Vein = T_N(NumDomRows+1+VesselRow:end);
    
    Tt(~DomTot) = NaN;
    Tb(~DomTot) = NaN;
    % Setting all values outside of the domain to NaN. This facilitates
    % visualisation of data with 'planecut'.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.