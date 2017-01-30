% this script calculates the AOD of each AeroCom model bin when mapped onto the bins of the four global models that calculated the radiative effect efficiency. This script also calculates the DRE from the size-resolved dust loading of each AeroCom model.

load('REE.mat'); %reading in the radiative effect efficiencies from the different models, written in write_REE
load('Vars_for_AeroCom_to_DRE_model_AOD.mat'); %the variables needed to calculate the DAOD from AeroCom models, mapped onto the DRE models, calculated in main_bootstrap_global_dust_cycle
load('Aerocom_data.mat'); %loading the AeroCom results, written in read_AeroCom_models
A_Earth = 510072000*10^6; %Earth surface area in m, from wiki/Earth
no_bootstrap_iters = size(Qe,1);
AeroCom_statistics = 'values'; %denotes whether the AeroCom statistics are calculated using a probability distribution or from the individual values. Options are 'values' and 'probdf'

%defining matrices
AeroCom_to_DRE_model_AOD_sp = zeros(no_AeroCom_models,max(no_AeroCom_bins),no_bootstrap_iters);
AeroCom_to_DRE_model_AOD_non_sp = zeros(no_AeroCom_models,max(no_AeroCom_bins),no_bootstrap_iters);
DAOD_AeroCom_sph = zeros(no_AeroCom_models,max(no_DRE_models),no_bootstrap_iters);
DAOD_AeroCom_m_DRE_k_sph = zeros(no_AeroCom_models,no_DRE_models,max(no_DRE_bins),no_bootstrap_iters);
DRE_AeroCom_m_DRE_k_sph_TOA_SW = zeros(no_AeroCom_models,no_DRE_models,max(no_DRE_bins),no_bootstrap_iters);
DRE_AeroCom_m_DRE_k_sph_TOA_LW = zeros(no_AeroCom_models,no_DRE_models,max(no_DRE_bins),no_bootstrap_iters);
DRE_AeroCom_m_DRE_k_sph_TOA_Tot = zeros(no_AeroCom_models,no_DRE_models,max(no_DRE_bins),no_bootstrap_iters);

%finding the i that corresponds to the AeroCom bin limits
for m=1:no_AeroCom_models %cycling over the seven AeroCom models
    for p=1:no_AeroCom_bins(m) %cycling over each model's particle bins
        if (D_AeroCom_lower(m,p) > min(D))
            D_AeroCom_i_bin_start(m,p) = find(D_AeroCom_lower(m,p)>=D+1e-8,1,'last')+1; %the 1e-8 is there to make the match at exactly 1 um for bin 2
        else
            D_AeroCom_i_bin_start(m,p) = 1;
        end
        if (D_AeroCom_upper(m,p) > min(D))
            D_AeroCom_i_bin_end(m,p) = find(D_AeroCom_upper(m,p)<=D,1,'first'); %find(D_AeroCom_upper(m,p)>=D,1,'last')+1;
        else
            D_AeroCom_i_bin_end(m,p) = 1;
        end    
    end %for, cycling over each model's particle bins
end %for, cycling over the seven AeroCom models

%finding the i that corresponds to the DRE model bin limits
for k=1:no_DRE_models %cycling over the seven DRE models
    for p=1:no_DRE_bins(k) %cycling over each model's particle bins
        if (D_DRE_lower(k,p) > min(D))
            D_DRE_i_bin_start(k,p) = find(D_DRE_lower(k,p)>=D+1e-8,1,'last')+1; %the 1e-8 is there to make the match at exactly 1 um for bin 2
        else
            D_DRE_i_bin_start(k,p) = 1;
        end
        if (D_DRE_upper(k,p) > min(D))
            D_DRE_i_bin_end(k,p) = find(D_DRE_upper(k,p)<=D,1,'first'); %find(D_DRE_upper(k,p)>=D,1,'last')+1;
        else
            D_DRE_i_bin_end(k,p) = 1;
        end    
    end %for, cycling over each model's particle bins
end %for, cycling over the DRE models

for m=1:no_AeroCom_models %cycling over the seven AeroCom models
    for i=1:size(D,2)
        %Calculating dM/dD for each model. cases:
        %1: below first bin or above last bin (not accounted for)
        %2: below midpoint of first bin: fit using first and second bins
        %3: above midpoint of last bin: fit using last and second-to-last bins
        %4: everewhere else
        if (D(i)<D_AeroCom_lower(m,1) || D(i)>D_AeroCom_upper(m,no_AeroCom_bins(m))+1e-10) %outside of the range of the AeroCom model
            dMdD_AeroCom_fct(m,i) = 0;
        elseif (D(i)<D_AeroCom_mean(m,2)) %below the midpoint of the second bin: fit using first and second bins
            k_start = 1; %starting bin for power law fit
            k_end = 2; %ending bin for power law fit
            x_fit = log(D_AeroCom_mean(m,k_start:k_end)); %x for power law fit of dM/dD to D
            y_fit = log(dMdD_AeroCom_bin(m,k_start:k_end)); %y for power law fit of dM/dD to D
            [intercept, slope] = linearfit(x_fit, y_fit, 0); %obtaining the power law's fitting coefficients
            dMdD_AeroCom_fct(m,i) = exp(intercept)*D(i)^slope; %the power law fit to dMdD (kg/m)
        elseif (D(i)>D_AeroCom_mean(m,no_AeroCom_bins(m)-1)) %above midpoint of second-to-last bin: fit using last and second-to-last bins
            k_start = no_AeroCom_bins(m)-1; %starting bin for power law fit
            k_end = no_AeroCom_bins(m); %ending bin for power law fit
            x_fit = log(D_AeroCom_mean(m,k_start:k_end)); %x for power law fit of dM/dD to D
            y_fit = log(dMdD_AeroCom_bin(m,k_start:k_end)); %y for power law fit of dM/dD to D
            [intercept, slope] = linearfit(x_fit, y_fit, 0); %obtaining the power law's fitting coefficients
            dMdD_AeroCom_fct(m,i) = exp(intercept)*D(i)^slope; %the power law fit to dMdD (kg/m)
        else %all other cases
            k_start = max(find(D(i)>D_AeroCom_mean(m,1:no_AeroCom_bins(m))));
            k_end = k_start + 1;
            x_fit = log(D_AeroCom_mean(m,k_start:k_end)); %x for power law fit of dM/dD to D
            y_fit = log(dMdD_AeroCom_bin(m,k_start:k_end)); %y for power law fit of dM/dD to D
            [intercept, slope] = linearfit(x_fit, y_fit, 0); %obtaining the power law's fitting coefficients
            dMdD_AeroCom_fct(m,i) = exp(intercept)*D(i)^slope; %the power law fit to dMdD (kg/m)
        end
    end
    dMdlnD_AeroCom_fct_kg(m,1:size(D,2)) = D_m.*dMdD_AeroCom_fct(m,1:size(D,2)); %dM/dlnD in kg; using that dM/dlnD = dD/dlnD * dM/dD = D*dM/dD
    dMdlnD_AeroCom_fct_Tg(m,1:size(D,2)) = 1e-9*dMdlnD_AeroCom_fct_kg(m,1:size(D,2)); %dM/dlnD in Tg
    
    %calculating the DAOD for each AeroCom bin size
    for p=1:no_AeroCom_bins(m)
        for j=1:no_bootstrap_iters
        	%calculate DAOD
            DAOD_AeroCom_sph(m,p,j) = (3/(2*rho_d(j)*A_Earth))*trapz(D_m(D_AeroCom_i_bin_start(m,p):D_AeroCom_i_bin_end(m,p)),dMdD_AeroCom_fct(m,D_AeroCom_i_bin_start(m,p):D_AeroCom_i_bin_end(m,p)).*Qe_Mie(D_AeroCom_i_bin_start(m,p):D_AeroCom_i_bin_end(m,p))./D_m(D_AeroCom_i_bin_start(m,p):D_AeroCom_i_bin_end(m,p))); %cannot use 'integral' because Qe is not a function. Only trust this when accurate is turned on.
        end
    end
    
    %calculating the tau_d of the psd of Aerocom model m, mapped onto the bin p of DRE model j
    for k=1:no_DRE_models %cycling over the DRE model results
    	for p=1:no_DRE_bins(k) %cycling over the particle bins of each DRE model
            for j=1:no_bootstrap_iters
                %calculate DAOD and DRE per bin in DRE model (spherical dust)
                DAOD_AeroCom_m_DRE_k_sph(m,k,p,j) = (3/(2*rho_d(j)*A_Earth))*trapz(D_m(D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p)),dMdD_AeroCom_fct(m,D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p)).*Qe_Mie(D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p))./D_m(D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p))); %cannot use 'integral' because Qe is not a function. Only trust this when accurate is turned on.
            end
            DRE_AeroCom_m_DRE_k_sph_TOA_SW(m,k,p,:) = DAOD_AeroCom_m_DRE_k_sph(m,k,p,:)*SW_REE_TOA(k,p);
            DRE_AeroCom_m_DRE_k_sph_TOA_LW(m,k,p,:) = DAOD_AeroCom_m_DRE_k_sph(m,k,p,:)*LW_REE_TOA(k,p);
            DRE_AeroCom_m_DRE_k_sph_TOA_Tot(m,k,p,:) = DAOD_AeroCom_m_DRE_k_sph(m,k,p,:)*(SW_REE_TOA(k,p)+LW_REE_TOA(k,p));
        end %for cycling over the particle bins of each DRE model
        %calculting DAOD and DRE for the up to 20 um bin
        if (D_AeroCom_upper(m,no_AeroCom_bins(m))>D_DRE_upper(k,no_DRE_bins(k))) %checking that the range of the AeroCom model exceeds that of the DRE model, in which case that excess range is extended to the PM20 bin
            for j=1:no_bootstrap_iters
            	%calculate DAOD and DRE per bin in DRE model (spherical dust)
                DAOD_AeroCom_m_DRE_k_sph_PM20(m,k,j) = (3/(2*rho_d(j)*A_Earth))*trapz(D_m(D_DRE_i_bin_end(k,no_DRE_bins(k)):end),dMdD_AeroCom_fct(m,D_DRE_i_bin_end(k,no_DRE_bins(k)):end).*Qe_Mie(D_DRE_i_bin_end(k,no_DRE_bins(k)):end)./D_m(D_DRE_i_bin_end(k,no_DRE_bins(k)):end)); %cannot use 'integral' because Qe is not a function. Only trust this when accurate is turned on.
            end
        else
            DAOD_AeroCom_m_DRE_k_sph_PM20(m,k,1:no_bootstrap_iters) = zeros(1,1,no_bootstrap_iters);
        end
        DRE_AeroCom_m_DRE_k_sph_TOA_SW_PM20(m,k,:) = DAOD_AeroCom_m_DRE_k_sph_PM20(m,k,:)*SW_REE_TOA_PM20;
        DRE_AeroCom_m_DRE_k_sph_TOA_LW_PM20(m,k,:) = DAOD_AeroCom_m_DRE_k_sph_PM20(m,k,:)*LW_REE_TOA_PM20;
        DRE_AeroCom_m_DRE_k_sph_TOA_Tot_PM20(m,k,:) = DAOD_AeroCom_m_DRE_k_sph_PM20(m,k,:)*(SW_REE_TOA_PM20+LW_REE_TOA_PM20);        
    end %for, cycling over the DRE model results
    m
end %cycling over AeroCom models

%constructing the arrays of DRE for each particle bin (p) for each DRE model (k) from the seven AeroCom
%PSDs projected onto that particle bin
for k=1:no_DRE_models
    for p=1:no_DRE_bins(k)
        DRE_AeroCom_sph_TOA_SW_bin(k,p,:) = reshape(squeeze(DRE_AeroCom_m_DRE_k_sph_TOA_SW(:,k,p,:)),1,[]);
        DRE_AeroCom_sph_TOA_LW_bin(k,p,:) = reshape(squeeze(DRE_AeroCom_m_DRE_k_sph_TOA_LW(:,k,p,:)),1,[]);
        DRE_AeroCom_sph_TOA_Tot_bin(k,p,:) = reshape(squeeze(DRE_AeroCom_m_DRE_k_sph_TOA_Tot(:,k,p,:)),1,[]);
    end
end

%averages of DAOD for AeroCom models
DAOD_AeroCom_sph_bin_mean = mean(DAOD_AeroCom_sph,3); %DAOD in AeroCom model per bin, spherical
total_DAOD_AeroCom_sph = sum(DAOD_AeroCom_sph_bin_mean,2); %total DAOD in AeroCom model, spherical

%calculating the DAOD for each AeroCom model bin, rescaled by the global DAOD reported in Table 3 of Huneeus et al. (2011); to be used in Fig. 2c
fid=fopen('DAOD_AeroCom_bins.txt','wt');
for m=1:no_AeroCom_models
    DAOD_AeroCom_sph_bin_scaled(m,:) = DAOD_AeroCom_sph_bin_mean(m,:).*(Tot_DAOD_AeroCom(m)/total_DAOD_AeroCom_sph(m));
    fprintf(fid,'D_mean%1.0f DAOD%1.0f ',m,m);
end
%writing the AeroCom DAOD data to a file
for k=1:max(no_AeroCom_bins)
    fprintf(fid,'\n');
    for m=1:no_AeroCom_models
        if (k<=no_AeroCom_bins(m))
            dDAODdlnD_Aerocom_sph_bin_scaled(m,k) = DAOD_AeroCom_sph_bin_scaled(m,k)/log(D_AeroCom_upper(m,k)/D_AeroCom_lower(m,k));
            fprintf(fid,'%2.3f %1.4e ',D_AeroCom_mean(m,k),dDAODdlnD_Aerocom_sph_bin_scaled(m,k));
        else %in this case, the model has fewer than k bins
            fprintf(fid,'0 0 ');
        end
    end
end
figure(1); clf; loglog(D,dMdlnD_AeroCom_fct_Tg); axis([0.2,20,0.01,30]); %plotting up dM/dlnD; should compare well against figure in paper

%averages of DAOD and DRE for AeroCom models, *mapped onto DRE models*, spherical
DAOD_AeroCom_DRE_model_sph_bin_mean = mean(DAOD_AeroCom_m_DRE_k_sph,4);
DRE_TOA_SW_AeroCom_sph_bin_mean = mean(DRE_AeroCom_m_DRE_k_sph_TOA_SW,4); 
DRE_TOA_LW_AeroCom_sph_bin_mean = mean(DRE_AeroCom_m_DRE_k_sph_TOA_LW,4);
DRE_TOA_Tot_AeroCom_sph_bin_mean = mean(DRE_AeroCom_m_DRE_k_sph_TOA_Tot,4);
DRE_TOA_SW_AeroCom_sph_PM20_mean = mean(DRE_AeroCom_m_DRE_k_sph_TOA_SW_PM20,3); 
DRE_TOA_LW_AeroCom_sph_PM20_mean = mean(DRE_AeroCom_m_DRE_k_sph_TOA_LW_PM20,3);
DRE_TOA_Tot_AeroCom_sph_PM20_mean = mean(DRE_AeroCom_m_DRE_k_sph_TOA_Tot_PM20,3);
DRE_TOA_SW_AeroCom_sph = sum(DRE_TOA_SW_AeroCom_sph_bin_mean,3);
DRE_TOA_LW_AeroCom_sph = sum(DRE_TOA_LW_AeroCom_sph_bin_mean,3);
%calculating the pdf of total DRE by cycling over all permutations of AeroCom model, that AeroCom model's SW DRE as predicted from each of the DRE models, and that AeroCom model's LW DRE as predicted from each of the DRE models
for m=1:no_AeroCom_models 
    for k=1:no_DRE_models
        for q=1:no_DRE_models
            DRE_TOA_Tot_AeroCom_sph(m,k,q)=DRE_TOA_SW_AeroCom_sph(m,k)+DRE_TOA_LW_AeroCom_sph(m,q);
        end
    end
end

%saving the DRE calculated by mapping each AeroCom model onto each DRE bin
for k=1:no_DRE_models
    filename = strcat('DRE_AeroCom_bins_model',int2str(k),'.txt');
    fid=fopen(filename,'wt');
    fprintf(fid,'Dmean DREsw DREswLowErr DREswHiErr DRElw DRElwLowErr DRElwHiErr DREtot DRElowErr DREhiErr \n');
    for p=1:no_DRE_bins(k) %writing out the data; the mean of the seven AeroCom estimates, and the lower error (mean-min) and upper error (max-mean)
        fprintf(fid,'%2.3f %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',sqrt(D_DRE_lower(k,p)*D_DRE_upper(k,p)),median(DRE_TOA_SW_AeroCom_sph_bin_mean(:,k,p)),median(DRE_TOA_SW_AeroCom_sph_bin_mean(:,k,p))-min(DRE_TOA_SW_AeroCom_sph_bin_mean(:,k,p)),max(DRE_TOA_SW_AeroCom_sph_bin_mean(:,k,p))-median(DRE_TOA_SW_AeroCom_sph_bin_mean(:,k,p)),median(DRE_TOA_LW_AeroCom_sph_bin_mean(:,k,p)),median(DRE_TOA_LW_AeroCom_sph_bin_mean(:,k,p))-min(DRE_TOA_LW_AeroCom_sph_bin_mean(:,k,p)),max(DRE_TOA_LW_AeroCom_sph_bin_mean(:,k,p))-median(DRE_TOA_LW_AeroCom_sph_bin_mean(:,k,p)),median(DRE_TOA_Tot_AeroCom_sph_bin_mean(:,k,p)),median(DRE_TOA_Tot_AeroCom_sph_bin_mean(:,k,p))-min(DRE_TOA_Tot_AeroCom_sph_bin_mean(:,k,p)),max(DRE_TOA_Tot_AeroCom_sph_bin_mean(:,k,p))-median(DRE_TOA_Tot_AeroCom_sph_bin_mean(:,k,p)));
    end %for
    %saving the PM20 bin
    fprintf(fid,'%2.3f %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',sqrt(D_DRE_upper(k,no_DRE_bins(k))*20),median(DRE_TOA_SW_AeroCom_sph_PM20_mean(:,k)),median(DRE_TOA_SW_AeroCom_sph_PM20_mean(:,k))-min(DRE_TOA_SW_AeroCom_sph_PM20_mean(:,k)),max(DRE_TOA_SW_AeroCom_sph_PM20_mean(:,k))-median(DRE_TOA_SW_AeroCom_sph_PM20_mean(:,k)),median(DRE_TOA_LW_AeroCom_sph_PM20_mean(:,k)),median(DRE_TOA_LW_AeroCom_sph_PM20_mean(:,k))-min(DRE_TOA_LW_AeroCom_sph_PM20_mean(:,k)),max(DRE_TOA_LW_AeroCom_sph_PM20_mean(:,k))-median(DRE_TOA_LW_AeroCom_sph_PM20_mean(:,k)),median(DRE_TOA_Tot_AeroCom_sph_PM20_mean(:,k)),median(DRE_TOA_Tot_AeroCom_sph_PM20_mean(:,k))-min(DRE_TOA_Tot_AeroCom_sph_PM20_mean(:,k)),max(DRE_TOA_Tot_AeroCom_sph_PM20_mean(:,k))-median(DRE_TOA_Tot_AeroCom_sph_PM20_mean(:,k)));
    fclose(fid);    
end
save('DRE_AeroCom_bins_model_data.mat','no_AeroCom_models','D_DRE_lower','D_DRE_upper','DRE_TOA_SW_AeroCom_sph_bin_mean','DRE_TOA_LW_AeroCom_sph_bin_mean','DRE_TOA_Tot_AeroCom_sph_bin_mean','DRE_TOA_SW_AeroCom_sph_PM20_mean','DRE_TOA_LW_AeroCom_sph_PM20_mean','DRE_TOA_Tot_AeroCom_sph_PM20_mean'); %saving model k data to matlab file

for m=1:no_AeroCom_models %correcting the DRE by scaling the DAOD to its AeroCom value in Table 3 of Huneeus et al. (2011)
    REE_TOA_SW_AeroCom_sph_scaled_temp(m,:) = DRE_TOA_SW_AeroCom_sph(m,:)/total_DAOD_AeroCom_sph(m);
    REE_TOA_LW_AeroCom_sph_scaled_temp(m,:) = DRE_TOA_LW_AeroCom_sph(m,:)/total_DAOD_AeroCom_sph(m);
    REE_TOA_Tot_AeroCom_sph_scaled_temp(m,:,:) = DRE_TOA_Tot_AeroCom_sph(m,:,:)/total_DAOD_AeroCom_sph(m);
    DRE_TOA_SW_AeroCom_sph_scaled_temp(m,:) = DRE_TOA_SW_AeroCom_sph(m,:)*(Tot_DAOD_AeroCom(m)/total_DAOD_AeroCom_sph(m));
    DRE_TOA_LW_AeroCom_sph_scaled_temp(m,:) = DRE_TOA_LW_AeroCom_sph(m,:)*(Tot_DAOD_AeroCom(m)/total_DAOD_AeroCom_sph(m));
    DRE_TOA_Tot_AeroCom_sph_scaled_temp(m,:,:) = DRE_TOA_Tot_AeroCom_sph(m,:,:)*(Tot_DAOD_AeroCom(m)/total_DAOD_AeroCom_sph(m));
end

%shaping the individual AeroCom estimates into a single array
no_SW_LW_DRE_estimates = no_AeroCom_models*no_DRE_models;
no_Tot_DRE_estimates = no_AeroCom_models*no_DRE_models*no_DRE_models;
REE_TOA_SW_AeroCom_sph_scaled = reshape(REE_TOA_SW_AeroCom_sph_scaled_temp,1,no_SW_LW_DRE_estimates);
REE_TOA_LW_AeroCom_sph_scaled = reshape(REE_TOA_LW_AeroCom_sph_scaled_temp,1,no_SW_LW_DRE_estimates);
REE_TOA_Tot_AeroCom_sph_scaled = reshape(REE_TOA_Tot_AeroCom_sph_scaled_temp,1,no_Tot_DRE_estimates);
DRE_TOA_SW_AeroCom_sph_scaled = reshape(DRE_TOA_SW_AeroCom_sph_scaled_temp,1,no_SW_LW_DRE_estimates);
DRE_TOA_LW_AeroCom_sph_scaled = reshape(DRE_TOA_LW_AeroCom_sph_scaled_temp,1,no_SW_LW_DRE_estimates);
DRE_TOA_Tot_AeroCom_sph_scaled = reshape(DRE_TOA_Tot_AeroCom_sph_scaled_temp,1,no_Tot_DRE_estimates);

%calculating the AeroCom statistics from the individual values; taking CI by eliminating the extreme points on either end
%calculating TOA statistics
SW_LW_outliers = 1; %number of outliers to be removed from either end for SW and LW
Tot_outliers = SW_LW_outliers*no_DRE_models; %number of outliers to be removed from either end for total DRE; larger than the value for SW_LW by a factor of no_DRE_outliers
A = sort(REE_TOA_SW_AeroCom_sph_scaled);
AeroCom_REE_TOA_SW_median = median(A);
AeroCom_REE_TOA_SW_upper_CI = A(end-SW_LW_outliers);
AeroCom_REE_TOA_SW_lower_CI = A(1+SW_LW_outliers);
A = sort(REE_TOA_LW_AeroCom_sph_scaled);
AeroCom_REE_TOA_LW_median = median(A);
AeroCom_REE_TOA_LW_upper_CI = A(end-SW_LW_outliers);
AeroCom_REE_TOA_LW_lower_CI = A(1+SW_LW_outliers);
A = sort(REE_TOA_Tot_AeroCom_sph_scaled);
AeroCom_REE_TOA_Tot_median = median(A);
AeroCom_REE_TOA_Tot_upper_CI = A(end-Tot_outliers);
AeroCom_REE_TOA_Tot_lower_CI = A(1+Tot_outliers);
A = sort(DRE_TOA_SW_AeroCom_sph_scaled);
AeroCom_DRE_TOA_SW_median = median(A);
AeroCom_DRE_TOA_SW_upper_CI = A(end-SW_LW_outliers);
AeroCom_DRE_TOA_SW_lower_CI = A(1+SW_LW_outliers);
A = sort(DRE_TOA_LW_AeroCom_sph_scaled);
AeroCom_DRE_TOA_LW_median = median(A);
AeroCom_DRE_TOA_LW_upper_CI = A(end-SW_LW_outliers);
AeroCom_DRE_TOA_LW_lower_CI = A(1+SW_LW_outliers);
A = sort(DRE_TOA_Tot_AeroCom_sph_scaled);
AeroCom_DRE_TOA_Tot_median = median(A);
AeroCom_DRE_TOA_Tot_upper_CI = A(end-Tot_outliers);
AeroCom_DRE_TOA_Tot_lower_CI = A(1+Tot_outliers);

save('AeroCom_DRE_results.mat','REE_TOA_SW_AeroCom_sph_scaled','REE_TOA_LW_AeroCom_sph_scaled','REE_TOA_Tot_AeroCom_sph_scaled','DRE_TOA_SW_AeroCom_sph_scaled','DRE_TOA_LW_AeroCom_sph_scaled','DRE_TOA_Tot_AeroCom_sph_scaled','AeroCom_REE_TOA_Tot_median','AeroCom_REE_TOA_Tot_upper_CI','AeroCom_REE_TOA_Tot_lower_CI','AeroCom_REE_TOA_SW_median','AeroCom_REE_TOA_SW_upper_CI','AeroCom_REE_TOA_SW_lower_CI','AeroCom_REE_TOA_LW_median','AeroCom_REE_TOA_LW_upper_CI','AeroCom_REE_TOA_LW_lower_CI','AeroCom_DRE_TOA_Tot_median','AeroCom_DRE_TOA_Tot_upper_CI','AeroCom_DRE_TOA_Tot_lower_CI','AeroCom_DRE_TOA_SW_median','AeroCom_DRE_TOA_SW_upper_CI','AeroCom_DRE_TOA_SW_lower_CI','AeroCom_DRE_TOA_LW_median','AeroCom_DRE_TOA_LW_upper_CI','AeroCom_DRE_TOA_LW_lower_CI','DRE_AeroCom_sph_TOA_SW_bin','DRE_AeroCom_sph_TOA_LW_bin','DRE_AeroCom_sph_TOA_Tot_bin');