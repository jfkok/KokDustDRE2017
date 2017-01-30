function main_calc_DRE (j_end, rerun_data)

%function that obtains constraints on the global dust DRE, using the size-resolved dust loading obtained through main_bootstrap_global_dust_cycle_size.m.

if (rerun_data==true)
    main_bootstrap_global_dust_cycle_size(j_end);
    calc_AeroCom_to_DRE_model_AOD;
end

%reading in dust AOD per particle size bin:
load('Model_AOD_per_bin.mat'); %reading in the DAOD constrained by analysis for the size range given by each model bin, calculated in main_bootstrap_global_dust_cycle_size
load('REE.mat'); %reading in the radiative effect efficiencies from the different models
load('AeroCom_DRE_results.mat'); %reading in the DRE for AeroCom models, written in calc_AeroCom_to_DRE_model_AOD
load('DRE_AeroCom_bins_model_data.mat'); %reads in the DRE per particle bin for the AeroCom models projected onto the DRE model particle bins, written in calc_AeroCom_to_DRE_model_AOD

%defining constants and REE for PM20 bin
A_Earth = 510072000*10^6; %Earth surface area in m, from wiki/Earth
DRE_pdf_bin_size = 0.02; %spacing of bins for DRE pdfs; in W/m2
%calculating TOA DRE for PM20 bin
for j = 1:j_end %randomly drawing values from the normal distribution describing the REE for the PM20 bin
	SW_REE_TOA_PM20_bootstrap(j)=Gaussian(SW_REE_TOA_PM20,SW_REE_TOA_PM20_SE);
    LW_REE_TOA_PM20_bootstrap(j)=Gaussian(LW_REE_TOA_PM20,LW_REE_TOA_PM20_SE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculating DRE for the CESM/GISS/GEOS-Chem simulations %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:no_DRE_models
    switch k
        case 1 %CESM
            fid_TOA = fopen('CESM_bin_TOA_DRE.txt','wt');
            fid_TOA_diff = fopen('CESM_bin_TOA_DRE_diff.txt','wt');
        case 2 %NASA GISS
            fid_TOA = fopen('GISS_bin_TOA_DRE.txt','wt');
            fid_TOA_diff = fopen('GISS_bin_TOA_DRE_diff.txt','wt');
        case 3 %GEOS-Chem
            fid_TOA = fopen('GEOS_Chem_bin_TOA_DRE.txt','wt');
            fid_TOA_diff = fopen('GEOS_Chem_bin_TOA_DRE_diff.txt','wt');
        case 4 %WRF-Chem
            fid_TOA = fopen('WRF_Chem_bin_TOA_DRE.txt','wt');
            fid_TOA_diff = fopen('WRF_Chem_bin_TOA_DRE_diff.txt','wt');
    end
    fprintf(fid_TOA,'Dlower  Dupper  SW_DRE_TOA_median  SW_DRE_TOA_lower_err  SW_DRE_TOA_upper_err  LW_DRE_TOA_median  LW_DRE_TOA_lower_err  LW_DRE_TOA_upper_err  Tot_DRE_TOA_median  Tot_DRE_TOA_lower_err  Tot_DRE_TOA_upper_err  \n');
    fprintf(fid_TOA_diff,'Dlower  Dupper  SW_DRE_TOA_median  SW_DRE_TOA_lower_err  SW_DRE_TOA_upper_err  LW_DRE_TOA_median  LW_DRE_TOA_lower_err  LW_DRE_TOA_upper_err  Tot_DRE_TOA_median  Tot_DRE_TOA_lower_err  Tot_DRE_TOA_upper_err  \n');
    for i=1:no_DRE_bins(k) %for figure S4
        %TOA, SW, our tau_d constraints
        TOA_SW_DRE_bin(k,i,:) = DRE_bin_AOD(k,:,i)*SW_REE_TOA(k,i);
        A=squeeze(sort(TOA_SW_DRE_bin(k,i,1:j_end)));
        TOA_SW_DRE_bin_median(k,i) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
        TOA_SW_DRE_bin_lower_CI(k,i) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
        TOA_SW_DRE_bin_upper_CI(k,i) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
        %TOA, SW, difference between AeroCom tau_d and our tau_d constraints
        TOA_SW_DRE_bin_diff = zeros(j_end*no_AeroCom_models,1);
        for q=1:no_AeroCom_models
            TOA_SW_DRE_bin_diff((1+(q-1)*j_end):q*j_end) = A-DRE_TOA_SW_AeroCom_sph_bin_mean(q,k,i);
        end
        A=squeeze(sort(TOA_SW_DRE_bin_diff));
        TOA_SW_DRE_diff_bin_median(k,i) = A(round(0.5*j_end*no_AeroCom_models)); %lower uncertainty value at -1 sigma
        TOA_SW_DRE_diff_bin_lower_CI(k,i) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end*no_AeroCom_models)))); %lower uncertainty value at -2 sigma
        TOA_SW_DRE_diff_bin_upper_CI(k,i) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end*no_AeroCom_models)))); %higher uncertainty value at +2 sigma        
        %TOA, LW, our tau_d constraints
        TOA_LW_DRE_bin(k,i,:) = DRE_bin_AOD(k,:,i)*LW_REE_TOA(k,i);        
        A=squeeze(sort(TOA_LW_DRE_bin(k,i,1:j_end)));
        TOA_LW_DRE_bin_median(k,i) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
        TOA_LW_DRE_bin_lower_CI(k,i) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
        TOA_LW_DRE_bin_upper_CI(k,i) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
        %TOA, LW, difference between AeroCom tau_d and our tau_d constraints
        TOA_LW_DRE_bin_diff = zeros(j_end*no_AeroCom_models,1);
        for q=1:no_AeroCom_models
            TOA_LW_DRE_bin_diff((1+(q-1)*j_end):q*j_end) = A-DRE_TOA_LW_AeroCom_sph_bin_mean(q,k,i);
        end
        A=squeeze(sort(TOA_LW_DRE_bin_diff));
        TOA_LW_DRE_diff_bin_median(k,i) = A(round(0.5*j_end*no_AeroCom_models)); %lower uncertainty value at -1 sigma
        TOA_LW_DRE_diff_bin_lower_CI(k,i) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end*no_AeroCom_models)))); %lower uncertainty value at -2 sigma
        TOA_LW_DRE_diff_bin_upper_CI(k,i) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end*no_AeroCom_models)))); %higher uncertainty value at +2 sigma        
        %TOA, Tot, our tau_d constraints
        TOA_Tot_DRE_bin(k,i,:) = DRE_bin_AOD(k,:,i)*(SW_REE_TOA(k,i)+LW_REE_TOA(k,i));
        A=squeeze(sort(TOA_Tot_DRE_bin(k,i,:)));
        TOA_Tot_DRE_bin_median(k,i) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
        TOA_Tot_DRE_bin_lower_CI(k,i) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
        TOA_Tot_DRE_bin_upper_CI(k,i) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
        %TOA, Tot, difference between AeroCom tau_d and our tau_d constraints
        TOA_Tot_DRE_bin_diff = zeros(j_end*no_AeroCom_models,1);
        for q=1:no_AeroCom_models
            TOA_Tot_DRE_bin_diff((1+(q-1)*j_end):q*j_end) = A-DRE_TOA_Tot_AeroCom_sph_bin_mean(q,k,i);
        end
        A=squeeze(sort(TOA_Tot_DRE_bin_diff));
        TOA_Tot_DRE_diff_bin_median(k,i) = A(round(0.5*j_end*no_AeroCom_models)); %lower uncertainty value at -1 sigma
        TOA_Tot_DRE_diff_bin_lower_CI(k,i) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end*no_AeroCom_models)))); %lower uncertainty value at -2 sigma
        TOA_Tot_DRE_diff_bin_upper_CI(k,i) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end*no_AeroCom_models)))); %higher uncertainty value at +2 sigma        
        %writing out the data for TOA
        fprintf(fid_TOA,'%2.3f  %2.3f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_DRE_lower(k,i), D_DRE_upper(k,i), TOA_SW_DRE_bin_median(k,i),TOA_SW_DRE_bin_median(k,i)-TOA_SW_DRE_bin_lower_CI(k,i),TOA_SW_DRE_bin_upper_CI(k,i)-TOA_SW_DRE_bin_median(k,i),TOA_LW_DRE_bin_median(k,i),TOA_LW_DRE_bin_median(k,i)-TOA_LW_DRE_bin_lower_CI(k,i),TOA_LW_DRE_bin_upper_CI(k,i)-TOA_LW_DRE_bin_median(k,i),TOA_Tot_DRE_bin_median(k,i),TOA_Tot_DRE_bin_median(k,i)-TOA_Tot_DRE_bin_lower_CI(k,i),TOA_Tot_DRE_bin_upper_CI(k,i)-TOA_Tot_DRE_bin_median(k,i));
        fprintf(fid_TOA_diff,'%2.3f  %2.3f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_DRE_lower(k,i), D_DRE_upper(k,i), TOA_SW_DRE_diff_bin_median(k,i),TOA_SW_DRE_diff_bin_median(k,i)-TOA_SW_DRE_diff_bin_lower_CI(k,i),TOA_SW_DRE_diff_bin_upper_CI(k,i)-TOA_SW_DRE_diff_bin_median(k,i),TOA_LW_DRE_diff_bin_median(k,i),TOA_LW_DRE_diff_bin_median(k,i)-TOA_LW_DRE_diff_bin_lower_CI(k,i),TOA_LW_DRE_diff_bin_upper_CI(k,i)-TOA_LW_DRE_diff_bin_median(k,i),TOA_Tot_DRE_diff_bin_median(k,i),TOA_Tot_DRE_diff_bin_median(k,i)-TOA_Tot_DRE_diff_bin_lower_CI(k,i),TOA_Tot_DRE_diff_bin_upper_CI(k,i)-TOA_Tot_DRE_diff_bin_median(k,i));
    end
    %PM20 bin, TOA, SW, our tau_d constraints
    TOA_SW_DRE_bin_PM20 = DRE_bin_AOD_PM20(k,:).*SW_REE_TOA_PM20_bootstrap;
    A=sort(TOA_SW_DRE_bin_PM20);
    TOA_SW_DRE_bin_PM20_median(k) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    TOA_SW_DRE_bin_PM20_lower_CI(k) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    TOA_SW_DRE_bin_PM20_upper_CI(k) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
    %PM20 bin, TOA, SW, difference between AeroCom tau_d and our tau_d constraints
    TOA_SW_DRE_bin_PM20_diff = zeros(j_end*no_AeroCom_models,1);
    for q=1:no_AeroCom_models
    	TOA_SW_DRE_bin_PM20_diff((1+(q-1)*j_end):q*j_end) = A-DRE_TOA_SW_AeroCom_sph_PM20_mean(q,k);
	end
    A=squeeze(sort(TOA_SW_DRE_bin_PM20_diff));
    TOA_SW_DRE_diff_bin_PM20_median(k) = A(round(0.5*j_end*no_AeroCom_models)); %lower uncertainty value at -1 sigma
    TOA_SW_DRE_diff_bin_PM20_lower_CI(k) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end*no_AeroCom_models)))); %lower uncertainty value at -2 sigma
	TOA_SW_DRE_diff_bin_PM20_upper_CI(k) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end*no_AeroCom_models)))); %higher uncertainty value at +2 sigma        
    %PM20 bin, TOA, LW, our tau_d constraints
    TOA_LW_DRE_bin_PM20 = DRE_bin_AOD_PM20(k,:).*LW_REE_TOA_PM20_bootstrap;
    A=sort(TOA_LW_DRE_bin_PM20);
    TOA_LW_DRE_bin_PM20_median(k) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    TOA_LW_DRE_bin_PM20_lower_CI(k) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    TOA_LW_DRE_bin_PM20_upper_CI(k) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
    %PM20 bin, TOA, LW, difference between AeroCom tau_d and our tau_d constraints
    TOA_LW_DRE_bin_PM20_diff = zeros(j_end*no_AeroCom_models,1);
    for q=1:no_AeroCom_models
    	TOA_LW_DRE_bin_PM20_diff((1+(q-1)*j_end):q*j_end) = A-DRE_TOA_LW_AeroCom_sph_PM20_mean(q,k);
	end
    A=squeeze(sort(TOA_LW_DRE_bin_PM20_diff));
    TOA_LW_DRE_diff_bin_PM20_median(k) = A(round(0.5*j_end*no_AeroCom_models)); %lower uncertainty value at -1 sigma
    TOA_LW_DRE_diff_bin_PM20_lower_CI(k) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end*no_AeroCom_models)))); %lower uncertainty value at -2 sigma
	TOA_LW_DRE_diff_bin_PM20_upper_CI(k) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end*no_AeroCom_models)))); %higher uncertainty value at +2 sigma        
    %PM20 bin, TOA, Tot, our tau_d constraints
    TOA_Tot_DRE_bin_PM20 = TOA_SW_DRE_bin_PM20+TOA_LW_DRE_bin_PM20;
    A=sort(TOA_Tot_DRE_bin_PM20);
    TOA_Tot_DRE_bin_PM20_median(k) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    TOA_Tot_DRE_bin_PM20_lower_CI(k) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    TOA_Tot_DRE_bin_PM20_upper_CI(k) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
    %PM20 bin, TOA, LW, difference between AeroCom tau_d and our tau_d constraints
    TOA_Tot_DRE_diff_bin_PM20 = TOA_SW_DRE_bin_PM20_diff+TOA_LW_DRE_bin_PM20_diff;
    A=squeeze(sort(TOA_Tot_DRE_diff_bin_PM20));
    TOA_Tot_DRE_diff_bin_PM20_median(k) = A(round(0.5*j_end*no_AeroCom_models)); %lower uncertainty value at -1 sigma
    TOA_Tot_DRE_diff_bin_PM20_lower_CI(k) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end*no_AeroCom_models)))); %lower uncertainty value at -2 sigma
	TOA_Tot_DRE_diff_bin_PM20_upper_CI(k) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end*no_AeroCom_models)))); %higher uncertainty value at +2 sigma        
    %writing out the TOA data
    fprintf(fid_TOA,'%2.3f  %2.3f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_DRE_upper(k,no_DRE_bins(k)), 20, TOA_SW_DRE_bin_PM20_median(k),TOA_SW_DRE_bin_PM20_median(k)-TOA_SW_DRE_bin_PM20_lower_CI(k),TOA_SW_DRE_bin_PM20_upper_CI(k)-TOA_SW_DRE_bin_PM20_median(k),TOA_LW_DRE_bin_PM20_median(k),TOA_LW_DRE_bin_PM20_median(k)-TOA_LW_DRE_bin_PM20_lower_CI(k),TOA_LW_DRE_bin_PM20_upper_CI(k)-TOA_LW_DRE_bin_PM20_median(k),TOA_Tot_DRE_bin_PM20_median(k),TOA_Tot_DRE_bin_PM20_median(k)-TOA_Tot_DRE_bin_PM20_lower_CI(k),TOA_Tot_DRE_bin_PM20_upper_CI(k)-TOA_Tot_DRE_bin_PM20_median(k));
    fprintf(fid_TOA_diff,'%2.3f  %2.3f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_DRE_upper(k,no_DRE_bins(k)), 20, TOA_SW_DRE_diff_bin_PM20_median(k),TOA_SW_DRE_diff_bin_PM20_median(k)-TOA_SW_DRE_diff_bin_PM20_lower_CI(k),TOA_SW_DRE_diff_bin_PM20_upper_CI(k)-TOA_SW_DRE_diff_bin_PM20_median(k),TOA_LW_DRE_diff_bin_PM20_median(k),TOA_LW_DRE_diff_bin_PM20_median(k)-TOA_LW_DRE_diff_bin_PM20_lower_CI(k),TOA_LW_DRE_diff_bin_PM20_upper_CI(k)-TOA_LW_DRE_diff_bin_PM20_median(k),TOA_Tot_DRE_diff_bin_PM20_median(k),TOA_Tot_DRE_diff_bin_PM20_median(k)-TOA_Tot_DRE_diff_bin_PM20_lower_CI(k),TOA_Tot_DRE_diff_bin_PM20_upper_CI(k)-TOA_Tot_DRE_diff_bin_PM20_median(k));
    fclose(fid_TOA);
    fclose(fid_TOA_diff);
    TOA_SW_DRE(k,:) = squeeze(sum(TOA_SW_DRE_bin(k,:,:)))'+TOA_SW_DRE_bin_PM20;
    TOA_LW_DRE(k,:) = squeeze(sum(TOA_LW_DRE_bin(k,:,:)))'+TOA_LW_DRE_bin_PM20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculating mean and CI DRE and REE for the CESM/GISS/GEOS-Chem/WRF-chem simulations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating mean and CI for TOA DRE
i_end = no_DRE_models*j_end;
TOA_SW_DRE_array=reshape(TOA_SW_DRE,1,[]);
TOA_LW_DRE_array=reshape(TOA_LW_DRE,1,[]);
for j=1:j_end %total DRE is random combination of LW and SW
    a_SW_rand(j) = ceil(no_DRE_models*rand(1));
    a_LW_rand(j) = ceil(no_DRE_models*rand(1));
    TOA_SW_DRE_array_save(j) = TOA_SW_DRE(a_SW_rand(j),j); %these arrays are saved to be written out, so the DRE can be decomposed into its LW and SW contributions
    TOA_LW_DRE_array_save(j) = TOA_LW_DRE(a_LW_rand(j),j); %these arrays are saved to be written out, so the DRE can be decomposed into its LW and SW contributions
    TOA_DRE_array(j) = TOA_SW_DRE(a_SW_rand(j),j)+TOA_LW_DRE(a_LW_rand(j),j); %Total DRE by randomly selecting a SW DRE from one of the DRE models and combining it with the LW DRE from a randomly-selected DRE model
end
figure(1); clf; plot(TOA_SW_DRE(1,:),'bs'); hold on; plot(TOA_SW_DRE(2,:),'gs'); plot(TOA_SW_DRE(3,:),'rs'); plot(TOA_SW_DRE(4,:),'ms');
figure(2); clf; plot(TOA_LW_DRE(1,:),'bs'); hold on; plot(TOA_LW_DRE(2,:),'gs'); plot(TOA_LW_DRE(3,:),'rs'); plot(TOA_LW_DRE(4,:),'ms');
TOA_DRE_correlated = TOA_SW_DRE+TOA_LW_DRE;
figure(3); clf; subplot(1,2,1); plot(TOA_DRE_array,'s'); axis([1 j_end -0.75 0.25]); hold on; subplot(1,2,2); plot(TOA_DRE_correlated(1,:),'bs'); hold on; plot(TOA_DRE_correlated(2,:),'gs'); plot(TOA_DRE_correlated(3,:),'rs'); plot(TOA_DRE_correlated(4,:),'ms'); axis([1 j_end -0.75 0.25]); 
A = sort(TOA_SW_DRE_array);
TOA_SW_DRE_median = A(round(0.5*i_end)); %lower uncertainty value at -1 sigma
TOA_SW_DRE_lower_CI = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*i_end)))); %lower uncertainty value at -2 sigma
TOA_SW_DRE_upper_CI = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*i_end)))); %higher uncertainty value at +2 sigma
A = sort(TOA_LW_DRE_array);
TOA_LW_DRE_median = A(round(0.5*i_end)); %lower uncertainty value at -1 sigma
TOA_LW_DRE_lower_CI = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*i_end)))); %lower uncertainty value at -2 sigma
TOA_LW_DRE_upper_CI = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*i_end)))); %higher uncertainty value at +2 sigma
A = sort(TOA_DRE_array);
prob_pos_DRE = size(find(A>0),2)/j_end;
TOA_DRE_median = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
TOA_DRE_lower_CI = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
TOA_DRE_upper_CI = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma
%calculating mean and CI for TOA REE
tau_d_array = repmat(tau_d_bootstrap,1,no_DRE_models);
TOA_SW_REE_array = TOA_SW_DRE_array./tau_d_array;
TOA_LW_REE_array = TOA_LW_DRE_array./tau_d_array;
TOA_SW_REE_array_save = TOA_SW_DRE_array_save./tau_d_bootstrap; %these arrays are saved to be written out, so the REE can be decomposed into its LW and SW contributions
TOA_LW_REE_array_save = TOA_LW_DRE_array_save./tau_d_bootstrap; %these arrays are saved to be written out, so the REE can be decomposed into its LW and SW contributions
TOA_REE_array = TOA_DRE_array./tau_d_bootstrap;
A = sort(TOA_SW_REE_array);
TOA_SW_REE_median = A(round(0.5*i_end)); %lower uncertainty value at -1 sigma
TOA_SW_REE_lower_CI = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*i_end)))); %lower uncertainty value at -2 sigma
TOA_SW_REE_upper_CI = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*i_end)))); %higher uncertainty value at +2 sigma
A = sort(TOA_LW_REE_array);
TOA_LW_REE_median = A(round(0.5*i_end)); %lower uncertainty value at -1 sigma
TOA_LW_REE_lower_CI = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*i_end)))); %lower uncertainty value at -2 sigma
TOA_LW_REE_upper_CI = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*i_end)))); %higher uncertainty value at +2 sigma
A = sort(TOA_REE_array);
TOA_REE_median = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
TOA_REE_lower_CI = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
TOA_REE_upper_CI = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %higher uncertainty value at +2 sigma

%comparing SW, LW, Tot DRE from AeroCom and from DRE models
[mean(DRE_TOA_SW_AeroCom_sph_scaled),mean(DRE_TOA_LW_AeroCom_sph_scaled),mean(DRE_TOA_Tot_AeroCom_sph_scaled)]
[TOA_SW_DRE_median,TOA_LW_DRE_median,TOA_DRE_median]

%comparing SW, LW, Tot REE from AeroCom and from DRE models
[mean(REE_TOA_SW_AeroCom_sph_scaled),mean(REE_TOA_LW_AeroCom_sph_scaled),mean(REE_TOA_Tot_AeroCom_sph_scaled)]
[TOA_SW_REE_median,TOA_LW_REE_median,TOA_REE_median]

%comparing total, with statistics, DRE
[AeroCom_DRE_TOA_Tot_median,AeroCom_DRE_TOA_Tot_lower_CI,AeroCom_DRE_TOA_Tot_upper_CI]
[TOA_DRE_median,TOA_DRE_lower_CI,TOA_DRE_upper_CI]

%writing out the REE arrays
fid=fopen('REE_array.txt','wt');
fprintf(fid,'SW_TOA LW_TOA Tot_TOA \n');
for i=1:j_end
    fprintf(fid,'%1.4e %1.4e %1.4e \n',TOA_SW_REE_array_save(i),TOA_LW_REE_array_save(i),TOA_REE_array(i));
end
fclose(fid);

%writing out the DRE arrays, to be read in by Origin
for k=1:no_DRE_models
    [temp_array_x,temp_array_prob] = bin(TOA_SW_DRE(k,:), DRE_pdf_bin_size, '', 0, 0);
    TOA_SW_DRE_x(k,1:size(temp_array_x,2))= temp_array_x; TOA_SW_DRE_prob(k,1:size(temp_array_x,2)) = temp_array_prob;
    [temp_array_x,temp_array_prob] = bin(TOA_LW_DRE(k,:), DRE_pdf_bin_size, '', 0, 0);
    TOA_LW_DRE_x(k,1:size(temp_array_x,2))= temp_array_x; TOA_LW_DRE_prob(k,1:size(temp_array_x,2)) = temp_array_prob;
end
TOA_SW_DRE_x(TOA_SW_DRE_x==0)=NaN; TOA_LW_DRE_x(TOA_LW_DRE_x==0)=NaN; %to prevent from plotting zeroes
%calculating the pdf of the combined DRE models
[TOA_SW_DRE_x_avg,TOA_SW_DRE_prob_avg] = bin(reshape(TOA_SW_DRE,1,[]), DRE_pdf_bin_size, '', 0, 0);
[TOA_LW_DRE_x_avg,TOA_LW_DRE_prob_avg] = bin(reshape(TOA_LW_DRE,1,[]), DRE_pdf_bin_size, '', 0, 0);
[TOA_Tot_DRE_x,TOA_Tot_DRE_prob] = bin(TOA_DRE_array, DRE_pdf_bin_size, '', 0, 0);
fid_SW=fopen('DRE_array_TOA_SW.txt','wt'); fid_LW=fopen('DRE_array_TOA_LW.txt','wt'); fid_Tot=fopen('DRE_array_TOA_Tot.txt','wt');
for k=1:no_DRE_models %writing the headers
    fprintf(fid_SW,'SW_TOA%1.0f_x SW_TOA%1.0f_prob ',k,k);
    fprintf(fid_LW,'LW_TOA%1.0f_x LW_TOA%1.0f_prob ',k,k);
    fprintf(fid_Tot,'Tot_TOA%1.0f_x Tot_TOA%1.0f_prob ',k,k);
end
fprintf(fid_SW,'SW_TOA_avg_x SW_TOA_avg_prob \n'); fprintf(fid_LW,'LW_TOA_avg_x LW_TOA_avg_prob \n'); fprintf(fid_Tot,'Tot_TOA_avg_x Tot_TOA_avg_prob \n');
for i=1:size(TOA_SW_DRE_x_avg,2)
    if (i<=size(TOA_SW_DRE_x,2)) %writing out the individual models
        for k=1:no_DRE_models
            fprintf(fid_SW,'%1.4e %1.4e ',TOA_SW_DRE_x(k,i),TOA_SW_DRE_prob(k,i));
        end
        fprintf(fid_SW,'%1.4e %1.4e ',TOA_SW_DRE_x_avg(i),TOA_SW_DRE_prob_avg(i));
    else
        for k=1:no_DRE_models
            fprintf(fid_SW,'NaN 0 ');
        end        
        fprintf(fid_SW,'%1.4e %1.4e ',TOA_SW_DRE_x_avg(i),TOA_SW_DRE_prob_avg(i));
    end        
    fprintf(fid_SW,'\n');
end
for i=1:size(TOA_LW_DRE_x_avg,2)
    if (i<=size(TOA_LW_DRE_x,2)) %writing out the individual models
        for k=1:no_DRE_models
            fprintf(fid_LW,'%1.4e %1.4e ',TOA_LW_DRE_x(k,i),TOA_LW_DRE_prob(k,i));
        end
        fprintf(fid_LW,'%1.4e %1.4e ',TOA_LW_DRE_x_avg(i),TOA_LW_DRE_prob_avg(i));
    else
        for k=1:no_DRE_models
            fprintf(fid_LW,'NaN 0 ');
        end
        fprintf(fid_LW,'%1.4e %1.4e ',TOA_LW_DRE_x_avg(i),TOA_LW_DRE_prob_avg(i));
    end        
    fprintf(fid_LW,'\n');
end
for i=1:size(TOA_Tot_DRE_x,2)
    fprintf(fid_Tot,'%1.4e %1.4e \n',TOA_Tot_DRE_x(i),TOA_Tot_DRE_prob(i));    
end
save('DRE_prob.mat','TOA_Tot_DRE_x','TOA_Tot_DRE_prob','TOA_SW_DRE_x_avg','TOA_SW_DRE_prob_avg','TOA_LW_DRE_x_avg','TOA_LW_DRE_prob_avg');
fclose(fid_SW); fclose(fid_LW); fclose(fid_Tot);

%writing out DRE, REE, their LW and SW components, and their upper and lower error bounds
%REE and DRE: SW, LW, Tot
fid=fopen('REE_TOA_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'x1 SW_AeroCom_REE SW_AeroCom_REE_lower_err SW_AeroCom_REE_upper_err x2 SW_REE SW_REE_lower_err SW_REE_upper_err x3 LW_AeroCom_REE LW_AeroCom_REE_lower_err LW_AeroCom_REE_upper_err x4 LW_REE LW_REE_lower_err LW_REE_upper_err x5 AeroCom_REE AeroCom_REE_lower_err AeroCom_REE_upper_err x6 REE REE_lower_err REE_upper_err \n');
fprintf(fid,'1 %1.4e %1.4e %1.4e 2 %1.4e %1.4e %1.4e 4 %1.4e %1.4e %1.4e 5 %1.4e %1.4e %1.4e 7 %1.4e %1.4e %1.4e 8 %1.4e %1.4e %1.4e \n',AeroCom_REE_TOA_SW_median,AeroCom_REE_TOA_SW_median-AeroCom_REE_TOA_SW_lower_CI,AeroCom_REE_TOA_SW_upper_CI-AeroCom_REE_TOA_SW_median,TOA_SW_REE_median,TOA_SW_REE_median-TOA_SW_REE_lower_CI,TOA_SW_REE_upper_CI-TOA_SW_REE_median,AeroCom_REE_TOA_LW_median,AeroCom_REE_TOA_LW_median-AeroCom_REE_TOA_LW_lower_CI,AeroCom_REE_TOA_LW_upper_CI-AeroCom_REE_TOA_LW_median,TOA_LW_REE_median,TOA_LW_REE_median-TOA_LW_REE_lower_CI,TOA_LW_REE_upper_CI-TOA_LW_REE_median,AeroCom_REE_TOA_Tot_median,AeroCom_REE_TOA_Tot_median-AeroCom_REE_TOA_Tot_lower_CI,AeroCom_REE_TOA_Tot_upper_CI-AeroCom_REE_TOA_Tot_median,TOA_REE_median,TOA_REE_median-TOA_REE_lower_CI,TOA_REE_upper_CI-TOA_REE_median);
fclose(fid);
fid=fopen('DRE_TOA_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'x1 SW_AeroCom_DRE SW_AeroCom_DRE_lower_err SW_AeroCom_DRE_upper_err x2 SW_DRE SW_DRE_lower_err SW_DRE_upper_err x3 LW_AeroCom_DRE LW_AeroCom_DRE_lower_err LW_AeroCom_DRE_upper_err x4 LW_DRE LW_DRE_lower_err LW_DRE_upper_err x5 AeroCom_DRE AeroCom_DRE_lower_err AeroCom_DRE_upper_err x6 DRE DRE_lower_err DRE_upper_err \n');
fprintf(fid,'1 %1.4e %1.4e %1.4e 2 %1.4e %1.4e %1.4e 4 %1.4e %1.4e %1.4e 5 %1.4e %1.4e %1.4e 7 %1.4e %1.4e %1.4e 8 %1.4e %1.4e %1.4e \n',AeroCom_DRE_TOA_SW_median,AeroCom_DRE_TOA_SW_median-AeroCom_DRE_TOA_SW_lower_CI,AeroCom_DRE_TOA_SW_upper_CI-AeroCom_DRE_TOA_SW_median,TOA_SW_DRE_median,TOA_SW_DRE_median-TOA_SW_DRE_lower_CI,TOA_SW_DRE_upper_CI-TOA_SW_DRE_median,AeroCom_DRE_TOA_LW_median,AeroCom_DRE_TOA_LW_median-AeroCom_DRE_TOA_LW_lower_CI,AeroCom_DRE_TOA_LW_upper_CI-AeroCom_DRE_TOA_LW_median,TOA_LW_DRE_median,TOA_LW_DRE_median-TOA_LW_DRE_lower_CI,TOA_LW_DRE_upper_CI-TOA_LW_DRE_median,AeroCom_DRE_TOA_Tot_median,AeroCom_DRE_TOA_Tot_median-AeroCom_DRE_TOA_Tot_lower_CI,AeroCom_DRE_TOA_Tot_upper_CI-AeroCom_DRE_TOA_Tot_median,TOA_DRE_median,TOA_DRE_median-TOA_DRE_lower_CI,TOA_DRE_upper_CI-TOA_DRE_median);
fclose(fid);
%REE and DRE: Tot only
fid=fopen('TOA_Tot_REE_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'AeroCom_REE AeroCom_REE_lower_err AeroCom_REE_upper_err REE REE_lower_err REE_upper_err \n');
fprintf(fid,'%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',AeroCom_REE_TOA_Tot_median,AeroCom_REE_TOA_Tot_median-AeroCom_REE_TOA_Tot_lower_CI,AeroCom_REE_TOA_Tot_upper_CI-AeroCom_REE_TOA_Tot_median,TOA_REE_median,TOA_REE_median-TOA_REE_lower_CI,TOA_REE_upper_CI-TOA_REE_median);
fclose(fid);
fid=fopen('TOA_Tot_DRE_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'AeroCom_DRE AeroCom_DRE_lower_err AeroCom_DRE_upper_err DRE DRE_lower_err DRE_upper_err \n');
fprintf(fid,'%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',AeroCom_DRE_TOA_Tot_median,AeroCom_DRE_TOA_Tot_median-AeroCom_DRE_TOA_Tot_lower_CI,AeroCom_DRE_TOA_Tot_upper_CI-AeroCom_DRE_TOA_Tot_median,TOA_DRE_median,TOA_DRE_median-TOA_DRE_lower_CI,TOA_DRE_upper_CI-TOA_DRE_median);
fclose(fid);
%REE and DRE: SW only:
fid=fopen('TOA_SW_REE_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'SW_AeroCom_REE SW_AeroCom_REE_lower_err SW_AeroCom_REE_upper_err SW_REE SW_REE_lower_err SW_REE_upper_err \n');
fprintf(fid,'%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',AeroCom_REE_TOA_SW_median,AeroCom_REE_TOA_SW_median-AeroCom_REE_TOA_SW_lower_CI,AeroCom_REE_TOA_SW_upper_CI-AeroCom_REE_TOA_SW_median,TOA_SW_REE_median,TOA_SW_REE_median-TOA_SW_REE_lower_CI,TOA_SW_REE_upper_CI-TOA_SW_REE_median);
fclose(fid);
fid=fopen('TOA_SW_DRE_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'SW_AeroCom_DRE SW_AeroCom_DRE_lower_err SW_AeroCom_DRE_upper_err SW_DRE SW_DRE_lower_err SW_DRE_upper_err \n');
fprintf(fid,'%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',AeroCom_DRE_TOA_SW_median,AeroCom_DRE_TOA_SW_median-AeroCom_DRE_TOA_SW_lower_CI,AeroCom_DRE_TOA_SW_upper_CI-AeroCom_DRE_TOA_SW_median,TOA_SW_DRE_median,TOA_SW_DRE_median-TOA_SW_DRE_lower_CI,TOA_SW_DRE_upper_CI-TOA_SW_DRE_median);
fclose(fid);
%REE and DRE: LW only
fid=fopen('TOA_LW_REE_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'LW_AeroCom_REE LW_AeroCom_REE_lower_err LW_AeroCom_REE_upper_err LW_REE LW_REE_lower_err LW_REE_upper_err \n');
fprintf(fid,'%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',AeroCom_REE_TOA_LW_median,AeroCom_REE_TOA_LW_median-AeroCom_REE_TOA_LW_lower_CI,AeroCom_REE_TOA_LW_upper_CI-AeroCom_REE_TOA_LW_median,TOA_LW_REE_median,TOA_LW_REE_median-TOA_LW_REE_lower_CI,TOA_LW_REE_upper_CI-TOA_LW_REE_median);
fclose(fid);
fid=fopen('TOA_LW_DRE_AeroCom_and_DRE_models.txt','wt');
fprintf(fid,'LW_AeroCom_DRE LW_AeroCom_DRE_lower_err LW_AeroCom_DRE_upper_err LW_DRE LW_DRE_lower_err LW_DRE_upper_err \n');
fprintf(fid,'%1.4e %1.4e %1.4e %1.4e %1.4e %1.4e \n',AeroCom_DRE_TOA_LW_median,AeroCom_DRE_TOA_LW_median-AeroCom_DRE_TOA_LW_lower_CI,AeroCom_DRE_TOA_LW_upper_CI-AeroCom_DRE_TOA_LW_median,TOA_LW_DRE_median,TOA_LW_DRE_median-TOA_LW_DRE_lower_CI,TOA_LW_DRE_upper_CI-TOA_LW_DRE_median);
fclose(fid);