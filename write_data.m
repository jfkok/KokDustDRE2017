%This script writes the important results from main_bootstrap_global_dust_cycle_size.m to files

%writing out lifetime bootstrap
fid = fopen('Tuncertainty.txt','wt');
fprintf(fid,'D  Max_likelihood Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
A=sort(T0.*exp(-0./D_dep));
T_max_likelihood = T0_max_likelihood*exp(-[0,D]/D_dep_max_likelihood);
T_avg(1) = A(ceil(0.5*j_end)); %upper uncertainty value at +1 sigma
T_neg_2sigma(1) = A(ceil(0.5*(1+erf(-2/sqrt(2)))*j_end)); %lower uncertainty value at -2 sigma
T_neg_1sigma(1) = A(ceil(0.5*(1+erf(-1/sqrt(2)))*j_end)); %lower uncertainty value at -1 sigma
T_pos_1sigma(1) = A(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)); %upper uncertainty value at +1 sigma
T_pos_2sigma(1) = A(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)); %upper uncertainty value at +2 sigma
fprintf(fid,'0.00  %1.4e  %1.4e  %1.4e  %1.4e %1.4e  %1.4e \n', T_max_likelihood(1), T_avg(1), T_neg_2sigma(1), T_neg_1sigma(1), T_pos_1sigma(1), T_pos_2sigma(1));
for p = 1:size(D,2)
    A=sort(T0.*exp(-D(p)./D_dep));
    T_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    T_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    T_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    T_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    T_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), T_max_likelihood(p+1), T_avg(p+1), T_neg_2sigma(p+1), T_neg_1sigma(p+1), T_pos_1sigma(p+1), T_pos_2sigma(p+1));
end
fclose(fid);

%writing out size-resolved emission flux bootstrap
fid = fopen('PSD_emit.txt','wt');
fprintf(fid,'D  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(PSD_lnV_emit(:,p));
    PSD_emit_avg(p+1) = A(round(0.5*j_end)); %mean bootstrap PSD
    PSD_emit_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_emit_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_emit_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_emit_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_emit_avg(p+1), PSD_emit_neg_2sigma(p+1), PSD_emit_neg_1sigma(p+1), PSD_emit_pos_1sigma(p+1), PSD_emit_pos_2sigma(p+1));
end
fclose(fid);

%writing out size-resolved normalized emission flux bootstrap
fid = fopen('PSD_emit_norm.txt','wt');
fprintf(fid,'D  Max_likelihood Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
PSD_lnV_max_likelihood = ([0,D]./Cv_max_likelihood).*(1+erf(log([0,D]./Ds_max_likelihood)/(sqrt(2)*log(sigma_s_max_likelihood)))).*exp(-([0,D]./lambda_max_likelihood).^3); %the PSD volume function from Kok, PNAS, 2011
for p = 1:size(D,2)
    A=sort(PSD_lnV_emit_norm(:,p));
    PSD_emit_avg_norm(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_emit_neg_2sigma_norm(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_emit_neg_1sigma_norm(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_emit_pos_1sigma_norm(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_emit_pos_2sigma_norm(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_lnV_max_likelihood(p+1), PSD_emit_avg_norm(p+1), PSD_emit_neg_2sigma_norm(p+1), PSD_emit_neg_1sigma_norm(p+1), PSD_emit_pos_1sigma_norm(p+1), PSD_emit_pos_2sigma_norm(p+1));
end
fclose(fid);

%writing out bins for size-resolved emission flux bootstrap
fid = fopen('PSD_emit_bins.txt','wt');
fprintf(fid,'D_min D_max  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D_min,2)
    A=sort(bin_emit(:,p));
    PSD_bin_emit_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_bin_emit_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_bin_emit_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_bin_emit_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_bin_emit_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_min(p), D_max(p), PSD_bin_emit_avg(p+1), PSD_bin_emit_neg_2sigma(p+1), PSD_bin_emit_neg_1sigma(p+1), PSD_bin_emit_pos_1sigma(p+1), PSD_bin_emit_pos_2sigma(p+1));
end
fclose(fid);

%writing out size-resolved dust volume loading bootstrap
fid = fopen('PSD_load.txt','wt');
fprintf(fid,'D  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(PSD_lnV_load(:,p));
    PSD_load_avg(p+1) = A(round(0.5*j_end)); %mean bootstrap PSD
    PSD_load_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_load_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_load_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_load_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_load_avg(p+1), PSD_load_neg_2sigma(p+1), PSD_load_neg_1sigma(p+1), PSD_load_pos_1sigma(p+1), PSD_load_pos_2sigma(p+1));
end
fclose(fid);

%writing out size-resolved normalized dust loading bootstrap
fid = fopen('PSD_load_norm.txt','wt');
func_int = @(D)(1./Cv_max_likelihood).*(1+erf(log(D./Ds_max_likelihood)/(sqrt(2)*log(sigma_s_max_likelihood)))).*exp(-(D./lambda_max_likelihood).^3)*T0_max_likelihood.*exp(-D./D_dep_max_likelihood); %function to integrate to get mean dust lifetime (see main paper)
normal_fact = integral(func_int,0.2,D_end); %the normalization factor
PSD_load_dVlnD_max_likelihood_norm = @(D)(D./Cv_max_likelihood).*(1+erf(log(D./Ds_max_likelihood)/(sqrt(2)*log(sigma_s_max_likelihood)))).*exp(-(D./lambda_max_likelihood).^3)*T0_max_likelihood.*exp(-D./D_dep_max_likelihood)/normal_fact; %function of the atmospheric size distribution
fprintf(fid,'D  Max_likelihood Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(PSD_lnV_load_norm(:,p));
    PSD_load_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_load_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_load_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_load_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_load_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_load_dVlnD_max_likelihood_norm(D(p)), PSD_load_avg(p+1), PSD_load_neg_2sigma(p+1), PSD_load_neg_1sigma(p+1), PSD_load_pos_1sigma(p+1), PSD_load_pos_2sigma(p+1));
end
fclose(fid);

%writing out bins for size-resolved dust loading bootstrap
fid = fopen('PSD_load_bins.txt','wt');
fprintf(fid,'D_min D_max  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D_min,2)
    A=sort(bin_load(:,p));
    PSD_bin_load_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_bin_load_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_bin_load_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_bin_load_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_bin_load_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_min(p), D_max(p), PSD_bin_load_avg(p+1), PSD_bin_load_neg_2sigma(p+1), PSD_bin_load_neg_1sigma(p+1), PSD_bin_load_pos_1sigma(p+1), PSD_bin_load_pos_2sigma(p+1));
end
fclose(fid);

%writing out size-resolved dust number loading bootstrap
fid = fopen('PSD_load_number_times1e-24.txt','wt');
fprintf(fid,'D  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(1e-24*PSD_lnN_load(:,p));
    PSD_load_N_avg(p+1) = A(round(0.5*j_end)); %mean bootstrap PSD
    PSD_load_N_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_load_N_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_load_N_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_load_N_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_load_N_avg(p+1), PSD_load_N_neg_2sigma(p+1), PSD_load_N_neg_1sigma(p+1), PSD_load_N_pos_1sigma(p+1), PSD_load_N_pos_2sigma(p+1));
end
fclose(fid);

%writing out size-resolved dust number per bin bootstrap
fid = fopen('PSD_bin_load_number_times1e-24.txt','wt');
fprintf(fid,'D_min D_max  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D_min,2)
    A=sort(1e-24*N_bin_load(:,p));
    PSD_N_bin_load_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_N_bin_load_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_N_bin_load_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_N_bin_load_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_N_bin_load_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_min(p), D_max(p), PSD_N_bin_load_avg(p+1), PSD_N_bin_load_neg_2sigma(p+1), PSD_N_bin_load_neg_1sigma(p+1), PSD_N_bin_load_pos_1sigma(p+1), PSD_N_bin_load_pos_2sigma(p+1));
end
fclose(fid);

%writing out the (not normalized) size-resolved global dust AOD bootstrap
fid = fopen('PSD_AOD.txt','wt');
fprintf(fid,'D  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(PSD_lnV_AOD(:,p));
    PSD_AOD_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_AOD_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_AOD_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_AOD_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_AOD_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_AOD_avg(p+1), PSD_AOD_neg_2sigma(p+1), PSD_AOD_neg_1sigma(p+1), PSD_AOD_pos_1sigma(p+1), PSD_AOD_pos_2sigma(p+1));
end
fclose(fid);

%writing out the 'normalized' (i.e., per 1 Tg of loading) size-resolved global dust AOD bootstrap
fid = fopen('PSD_AOD_norm.txt','wt');
fprintf(fid,'D  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(PSD_lnV_AOD_norm(:,p));
    PSD_AOD_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_AOD_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_AOD_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_AOD_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_AOD_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), PSD_AOD_avg(p+1), PSD_AOD_neg_2sigma(p+1), PSD_AOD_neg_1sigma(p+1), PSD_AOD_pos_1sigma(p+1), PSD_AOD_pos_2sigma(p+1));
end
fclose(fid);

%writing out bins for size-resolved global dust AOD bootstrap
fid = fopen('PSD_AOD_bins.txt','wt');
fprintf(fid,'D_min D_max  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma \n');
for p = 1:size(D_min,2)
    A=sort(bin_AOD(:,p));
    PSD_bin_AOD_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    PSD_bin_AOD_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    PSD_bin_AOD_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    PSD_bin_AOD_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    PSD_bin_AOD_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D_min(p), D_max(p), PSD_bin_AOD_avg(p+1), PSD_bin_AOD_neg_2sigma(p+1), PSD_bin_AOD_neg_1sigma(p+1), PSD_bin_AOD_pos_1sigma(p+1), PSD_bin_AOD_pos_2sigma(p+1));
end
fclose(fid);

%writing out extinction efficiency with bootstrapped uncertainty
fid = fopen('Qext.txt','wt');
fprintf(fid,'D Qext_mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma Qe_Mie \n');
for p = 1:size(D,2)
    A=sort(Qe(:,p));
    Qext_avg(p) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    Qext_neg_2sigma(p) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    Qext_neg_1sigma(p) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    Qext_pos_1sigma(p) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    Qext_pos_2sigma(p) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma    
    fprintf(fid,'%1.5f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), Qext_avg(p), Qext_neg_2sigma(p), Qext_neg_1sigma(p), Qext_pos_1sigma(p), Qext_pos_2sigma(p), Qe_Mie(p));
end
fclose(fid);

%writing out contributions to optical depth as a function of particle size
fid = fopen('DAOD_vs_D.txt','wt');
fprintf(fid,'D  Mean  Neg_2sigma  Neg_1sigma  Pos_1sigma Pos_2sigma  Cum_mean  Cum_neg_2sigma  Cum_neg_1sigma  Cum_pos_1sigma Cum_pos_2sigma \n');
for p = 1:size(D,2)
    A=sort(AOD_rel_contribution(:,p));
    DAOD_vs_D_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    DAOD_vs_D_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    DAOD_vs_D_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    DAOD_vs_D_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    DAOD_vs_D_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    A=sort(AOD_rel_contribution_cum(:,p));
    DAOD_vs_D_cum_avg(p+1) = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
    DAOD_vs_D_cum_neg_2sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
    DAOD_vs_D_cum_neg_1sigma(p+1) = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
    DAOD_vs_D_cum_pos_1sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    DAOD_vs_D_cum_pos_2sigma(p+1) = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
    fprintf(fid,'%1.4f  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e %1.4e  %1.4e  %1.4e  %1.4e  %1.4e \n', D(p), DAOD_vs_D_avg(p+1), DAOD_vs_D_neg_2sigma(p+1), DAOD_vs_D_neg_1sigma(p+1), DAOD_vs_D_pos_1sigma(p+1), DAOD_vs_D_pos_2sigma(p+1), DAOD_vs_D_cum_avg(p+1), DAOD_vs_D_cum_neg_2sigma(p+1), DAOD_vs_D_cum_neg_1sigma(p+1), DAOD_vs_D_cum_pos_1sigma(p+1), DAOD_vs_D_cum_pos_2sigma(p+1));
end
fclose(fid);