function main_bootstrap_global_dust_cycle_size (j_end)

%This function uses a bootstrap procedure to constrain the size-resolved global dust loading, using Eqs. (1) – (4) in Kok et al. (2017)

global D_end T0_range_warning D_dep_range_warning Ds_range_warning sigma_s_range_warning lambda_range_warning calculate_likelihood_function_gradient

%constants and DAOD
lambda550 = 0.550; %wavelength at which extinction is measured in um
A_Earth = 510072000*10^6; %Earth surface area in m, from wiki/Earth
daysInYears = 365;
load('DAOD.mat'); %loading the pdf of the DAOD from Ridley et al. (2016)
load('Aerocom_data.mat'); %loading the results from AeroCom models; j = 1 = CAM, 2 = GISS, 3 = GOCART, 4 = MATCH, 5 = MOZART, 6 = IMPACT, 7 = LOA
load('REE.mat'); %loading details on the DRE model results; j = 1 = CESM, 2 = GISS, 3 = GEOS-Chem, 4 = WRF-Chem

D_end = 20; %the upper limit on the size range
mean_rho_d = 2500; %density of dust particles, in kg/m3
error_rho_d = 200; %error in dust density
options = optimoptions(@fminunc,'GradObj','off');

%switches
recalculate_lifetime_and_PSD_fits = false; %sets whether to recalculate the fits to the PSD and lifetime data
calculate_likelihood_function_gradient = true; %sets whether to use the local gradient for the algorithm in the max likelihood method

%results from lifetime and PSD fits
if (recalculate_lifetime_and_PSD_fits)
    [T_intercept, T_slope, T_SE_intercept, T_SE_slope, T_rho, T0_max_likelihood, D_dep_max_likelihood, Cn, Ds, sigma_s, lambda, lambda_defined, PSD_covariance_matrix, Ds_max_likelihood, sigma_s_max_likelihood, lambda_max_likelihood] = fit_PSD_and_lifetime_functions;
end
load('lifetime_and_PSD_fits.mat'); %loads the fitting results from the data file
Ds_err = sqrt(squeeze(PSD_covariance_matrix(:,1,1))');
sigma_s_err = sqrt(squeeze(PSD_covariance_matrix(:,2,2))');
lambda_err = sqrt(squeeze(PSD_covariance_matrix(:,3,3))');
rho_PSD = PSD_covariance_matrix(:,1,2)'./(Ds_err.*sigma_s_err);
no_T_data_sets = size(T_intercept,2); no_PSD_data_sets = size(Ds_err,2);

%dust optical properties
load('TAMU_data.mat'); %reading in the extinction coefficient info from the TAMU database
load('TAMU_database_parameters.mat'); %contains the parameters for which the TAMU database was generated
n_mean = 1.53; %mean real refractive index
n_err = 0.03; %error in real refractive index
n_cum_prob = 0.5*(1+erf((n_database-n_mean)/(sqrt(2)*n_err)));
log_k_mean = -2.5; %imaginary refractive index; corresponds to -0.003i
log_k_err = 0.3; %error in imaginary refractive index; corresponds to 1 sd range of 0.0013 - 0.0079
log_k_cum_prob = 0.5*(1+erf((log10(k_database)-log_k_mean)/(sqrt(2)*log_k_err)));
AR_sigma_mean = 0.6; %the geometric standard deviation of the lognormal distribution that describes the aspect ratio
AR_sigma_error = 0.2; %error in the geometric standard deviation
AR_sigma_min = 0.1; %the minimum geometric standard deviation 
AR_median_mean = 1.7; %the median of the lognormal distribution that describes the aspect ratio
AR_median_error = 0.2; %the error in the median of the lognormal distribution that describes the aspect ratio
AR_median_min = 1.1; %the minimum median of the lognormal distribution that describes the aspect ratio
p_exp=1.6075; %exponent used in calculating approximate surface area of triaxial ellipsoid, from http://en.wikipedia.org/wiki/Ellipsoid#Surface_area
T0_range_warning = 0; D_dep_range_warning = 0; Ds_range_warning = 0; sigma_s_range_warning = 0; lambda_range_warning = 0; %initializing the global parameters keeping track of whether the parameter space sampled is too small
no_failed_max_likelihood_T = 0; %keeps track of number of times max likelihood method fails for lifetime;
no_failed_max_likelihood_PSD = 0; %keeps track of number of times max likelihood method fails for PSD;

%Obtaining the extinction coefficient from Mie theory
x = 3.14159*D/lambda550; %size parameter based on geometric diameter
for p=1:size(x,2)
    result = Mie(n_mean+10^log_k_mean*1i,x(p));       
    Qe_Mie(p) = result(1); %extinction coefficient from Mie theory
end
D_m = 1e-6*D; %particle size in m (D is in um)

%bins for bar charts
D_min = [0.1, 0.2, 0.5, 1, 2, 5, 10]; %geometric diameter in um
D_max = [0.2, 0.5, 1, 2, 5, 10, 20]; %geometric diameter in um
for p=1:size(D_min,2)
    if (p==1)
        D_i_bin_start(p) = 1;
    else
        D_i_bin_start(p) = find(D_min(p)>=D,1,'last')+1;
    end
    D_i_bin_end(p) = find(D_max(p)<=D,1,'first');
end

for k=1:no_DRE_models %cycling over the DRE models
    for p=1:no_DRE_bins(k) %finding the i that corresponds to the bin limits
        if (p==1)
            D_DRE_i_bin_start(k,p) = 1;
        else
            D_DRE_i_bin_start(k,p) = find(D_DRE_lower(k,p)>=D+1e-8,1,'last')+1; %the 1e-8 is there to make the match at exactly 1 um for bin 2
        end
        D_DRE_i_bin_end(k,p) = find(D_DRE_upper(k,p)<=D,1,'first');    
    end %for, finding the i that corresponds to the bin limits
end %for, cycling over the DRE models
%Size bins for GEOS-Chem

%initializing matrices and arrays
rho_d = zeros(1,j_end); AR_sigma = zeros(1,j_end); AR_median = zeros(1,j_end); AR_log_mu = zeros(1,j_end); A_ellipsoid_approx = zeros(1,j_end); surface_area_enhancement_over_spherical = zeros(1,j_end);
median_aspect_ratio = zeros(1,j_end); expected_median_aspect_ratio = zeros(1,j_end); surface_to_volume_ratio = zeros(1,j_end); geometric_ext_efficiency = zeros(1,j_end); 
lambda_bootstrap = zeros(1,j_end); lambda_bootstrap_alt = zeros(1,j_end); fval = zeros(1,j_end); exitflag = zeros(1,j_end); sigma_s_bootstrap = zeros(1,j_end); Ds_bootstrap = zeros(1,j_end);
A_min = zeros(1,j_end); Cv_bootstrap = zeros(1,j_end); Cn_bootstrap = zeros(1,j_end); T0 = zeros(1,j_end); D_dep = zeros(1,j_end); denominator = zeros(1,j_end);
denominator_D10 = zeros(1,j_end); F_eff = zeros(1,j_end); integral_Ltot = zeros(1,j_end); integral_Ltot_D10 = zeros(1,j_end); L_eff = zeros(1,j_end); tau_d_bootstrap = zeros(1,j_end);
F_emit = zeros(1,j_end); L_atm = zeros(1,j_end); column_burden = zeros(1,j_end); MEE = zeros(1,j_end); scattering_enhancement = zeros(1,j_end); N_IN = zeros(1,j_end);
clay_emit_fr = zeros(1,j_end); clay_load_fr = zeros(1,j_end); clay_N_load_fr = zeros(1,j_end); clay_AOD_fr = zeros(1,j_end); F_emit_D10 = zeros(1,j_end); L_atm_D10 = zeros(1,j_end);
column_burden_D10 = zeros(1,j_end); MEE_D10 = zeros(1,j_end); L_atm_PM1_geom = zeros(1,j_end); L_atm_PM25_geom = zeros(1,j_end); 

Qe = zeros(j_end,size(x,2)); Qe_TAMU = zeros(j_end,size(x,2)); g_TAMU = zeros(j_end,size(x,2)); sample_array_PSD = zeros(j_end,size(no_PSD_data_sets,2)); sample_array_T = zeros(j_end,size(no_T_data_sets,2)); X_min_PSD = zeros(j_end,2);
PSD_V = zeros(j_end,size(D,2)); PSD_N = zeros(j_end,size(D,2)); X_min_T = zeros(j_end,2); T = zeros(j_end,size(D,2)); PSD_lnV_emit = zeros(j_end,size(D,2)); PSD_lnV_emit_norm = zeros(j_end,size(D,2));
PSD_lnV_load = zeros(j_end,size(D,2)); PSD_lnV_load_norm = zeros(j_end,size(D,2)); PSD_lnN_load = zeros(j_end,size(D,2)); PSD_lnV_AOD = zeros(j_end,size(D,2)); 
PSD_lnV_AOD_norm = zeros(j_end,size(D,2)); AOD_rel_contribution = zeros(j_end,size(D,2)); AOD_rel_contribution_cum = zeros(j_end,size(D,2)); bin_emit = zeros(j_end,size(D,2)); 
bin_load = zeros(j_end,size(D,2)); N_bin_load = zeros(j_end,size(D,2)); m_part = zeros(j_end,size(D,2)); N_bin_load_alt = zeros(j_end,size(D,2)); bin_AOD = zeros(j_end,size(D,2)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%     Doing the resampling procedure     %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:j_end %looping over the resampling iterations
    if (mod(j,10)==0)
        j
    end
    %performing the dust density bootstrap
    rho_d(j) = Gaussian(mean_rho_d,error_rho_d);
    
    %performing the extinction efficiency resampling
    rand_no = rand(1); [A,i_min] = min(abs(n_cum_prob-rand_no)); n(j) = n_database(i_min); %drawing a realization from the normal distribution and then picking the value in the TAMU data set that is closest to it (in cum prob distr space)
    rand_no = rand(1); [A,i_min] = min(abs(log_k_cum_prob-rand_no)); k(j) = round(1e4*k_database(i_min))/1e4; %drawing a realization from the normal distribution and then picking the value in the TAMU data set that is closest to it (in cum prob distr space)
    AR_sigma(j) = max([AR_sigma_min,Gaussian(AR_sigma_mean,AR_sigma_error)]);
    AR_median(j) = max([AR_median_min,Gaussian(AR_median_mean,AR_median_error)]);
    AR_log_mu(j) = log(AR_median(j)-1);
    PSD_AR = @(AR) (1./(sqrt(2*pi)*(AR-1)*AR_sigma(j))).*exp(-0.5*((log(AR-1)-AR_log_mu(j))/AR_sigma(j)).^2); %PSD of aspect ratio, used to weigh the likelihood of occurrence of any particular aspect ratios into Qext. See https://en.wikipedia.org/wiki/Log-normal_distribution and Kandler et al. (2007, 2011) for exact equation
    for q=1:size(x,2)
        use_array = find(Vol_eq_size_parameter<1.01*x(q) & Vol_eq_size_parameter>0.99*x(q) & n_TAMU<1.0001*n(j) & n_TAMU>0.9999*n(j) & k_TAMU<1.0001*k(j) & k_TAMU>0.9999*k(j));
        Qe_TAMU(j,q) = sum((PSD_AR(aspect_ratio(use_array)).*Qext_wrt_VES_area(use_array)))/sum(PSD_AR(aspect_ratio(use_array))); %weighted mean extinction efficiency (collision cross section dividid by surface area of volume equivalent sphere)
    end
    
    %calculating the geometric extinction efficiency from the particle shape. Since the pdf of aspect ratio is not dependent on size, do not need to loop over particle sizes, and can do only once for each bootstrap sample
    a = (D(q)/2)*(aspect_ratio(use_array)./HLR).^(1/3); %semi-major axis in um
    b = a./aspect_ratio(use_array); %semi-minor axis in um
    c = a.*HLR; %half height in um
    A_ellipsoid_approx(j) = PSD_AR(aspect_ratio(use_array))'*4*3.14159*(((a.*b).^p_exp+(a.*c).^p_exp+(b.*c).^p_exp)/3).^(1/p_exp)/sum(PSD_AR(aspect_ratio(use_array))'); %from http://en.wikipedia.org/wiki/Ellipsoid#Surface_area; in um2
    surface_area_enhancement_over_spherical(j) = A_ellipsoid_approx(j)/(pi*D(q)^2); %enhancement in surface area over the spherical case
    median_aspect_ratio(j) = exp(sum(PSD_AR(aspect_ratio(use_array)).*log(aspect_ratio(use_array)))/sum(PSD_AR(aspect_ratio(use_array)))); %the mean in log space - should correspond to median
    expected_median_aspect_ratio(j) = exp(AR_log_mu(j)+0.5*AR_sigma(j)^2)+1; %expected median from the definition of the lognormal
    surface_to_volume_ratio(j) = A_ellipsoid_approx(j)/((pi/6)*D(q)^3); %in um-1
    geometric_ext_efficiency(j) = 2*A_ellipsoid_approx(j)./(pi*D(q)^2); %geometric extinction efficiency is two times the area, normalized by the area of a sphere with the same geometric diameter
    Qe(j,:) = Qe_TAMU(j,:);
    
    %performing the PSD resampling
    found_max_likelihood_PSD = false;
    while (~found_max_likelihood_PSD)
        try
            sample_array_PSD(j,1:no_PSD_data_sets) = randi([1 no_PSD_data_sets],1,no_PSD_data_sets);
            if (max(lambda_defined(sample_array_PSD(j,1:no_PSD_data_sets)))==0) %in this case, none of the selected data sets determined lambda
                lambda_bootstrap(j) = lambda_max_likelihood; 
            else
                A = lambda(sample_array_PSD(j,1:no_PSD_data_sets)); lambda_fit = A(lambda_defined(sample_array_PSD(j,1:no_PSD_data_sets))); %values of lambda_err to be used for determining its bootstrap value
                A = lambda_err(sample_array_PSD(j,1:no_PSD_data_sets)); SE_lambda_fit = A(lambda_defined(sample_array_PSD(j,1:no_PSD_data_sets))); %values of lambda_err to be used for determining its bootstrap value
                X0 = [mean(lambda_fit),std(lambda_fit)];
                L = @(X)1/(prod((1./sqrt(2*pi*(SE_lambda_fit.^2+X(2)^2))).*exp(-((lambda_fit-X(1)).^2)./(2*(SE_lambda_fit.^2+X(2)^2)))));
                [lambda_min(j,1:2),fval(j),exitflag(j)] = fminunc(L,X0,options); %finding the parameters that maximize the lifetime likelihood function
                lambda_bootstrap(j) = lambda_min(j,1); tau_lambda = lambda_min(j,2); %obtaining the fitting parameters
                lambda_bootstrap_alt(j) = sum(lambda_fit./(SE_lambda_fit.^2))./sum(1./(SE_lambda_fit.^2)); %calculating the weighted mean of the bootstrap lambda values                                
            end
            Ds_fit = Ds(sample_array_PSD(j,:)); sigma_s_fit = sigma_s(sample_array_PSD(j,:)); SE_Ds_fit = Ds_err(sample_array_PSD(j,:)); SE_sigma_s_fit = sigma_s_err(sample_array_PSD(j,:)); rho = rho_PSD(sample_array_PSD(j,:)); %setting the global fitting parameters equal to those for this particular bootstrap for the lifetime, in order to do the maximum likelihood estimation
            X0 =[mean(Ds_fit),mean(sigma_s_fit)]; %initial parameters vector
            %defining the variables for the likelihood function
            var_a = Ds_fit; SE_var_a = SE_Ds_fit;
            var_b = sigma_s_fit; SE_var_b = SE_sigma_s_fit;
            L = @(X)1/prod((1./(2*pi*sqrt((SE_var_a.^2+tau_Ds^2).*(SE_var_b.^2+tau_sigma_s^2).*(1-((rho.*SE_var_a.*SE_var_b+rho_syst_estimate*tau_Ds*tau_sigma_s)./(sqrt(SE_var_a.^2+tau_Ds^2).*sqrt(SE_var_b.^2+tau_sigma_s^2))).^2)))).*exp(-((var_a-X(1)).^2./(SE_var_a.^2+tau_Ds^2)+(var_b-X(2)).^2./(SE_var_b.^2+tau_sigma_s^2)-2.*((rho.*SE_var_a.*SE_var_b+rho_syst_estimate*tau_Ds*tau_sigma_s)./(sqrt(SE_var_a.^2+tau_Ds^2).*sqrt(SE_var_b.^2+tau_sigma_s^2))).*(var_a-X(1)).*(var_b-X(2))./(sqrt(SE_var_a.^2+tau_Ds^2).*sqrt(SE_var_b.^2+tau_sigma_s^2)))./(2.*(1-((rho.*SE_var_a.*SE_var_b+rho_syst_estimate*tau_Ds*tau_sigma_s)./(sqrt(SE_var_a.^2+tau_Ds^2).*sqrt(SE_var_b.^2+tau_sigma_s^2))).^2))));
            [X_min_PSD(j,1:2),fval(j),exitflag(j)] = fminunc(L,X0,options); %finding the parameters that maximize the lifetime likelihood function
            Ds_bootstrap(j) = X_min_PSD(j,1); sigma_s_bootstrap(j) = X_min_PSD(j,2);
        if (fval(j)>1e-25)
            found_max_likelihood_PSD = true;
        else
            error('problem in PSD maximum likelihood method: fval too low')            
        end
        catch
            no_failed_max_likelihood_PSD = no_failed_max_likelihood_PSD+1;
            'error in PSD maximum likelihood method'
            failed_sample_array_PSD(no_failed_max_likelihood_PSD,1:no_PSD_data_sets) = sample_array_PSD(j,1:no_PSD_data_sets);
        end
    end
    PSD = @(D) (1./D.^2).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %the PSD function from Kok, PNAS, 2011
    X0 = 0; %initial guess for A
    X_data = log10(D(D>2 & D<10));
    Y_data = log10(PSD(D(D>2 & D<10)));
    A_min(j) = lsqcurvefit(@dustPSD_power_law_fit, X0, X_data, Y_data); %proportionality factor of -2 power law with which function will be normalized
    y_int = @(D) (1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %the function to be integrated
    Cv_bootstrap(j) = integral(y_int,0.2,D_end); %calculating the normalization constant for the volume PSD
    y_int_N = @(D) (1./D.^3).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %the function to be integrated
    Cn_bootstrap(j) = integral(y_int_N,0.2,D_end); %calculating the normalization constant for the number PSD
    PSD_V(j,:) = (1./Cv_bootstrap(j)).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D./lambda_bootstrap(j)).^3); %the PSD volume function from Kok, PNAS, 2011
    PSD_N(j,:) = (1./(D.^3*Cn_bootstrap(j))).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D./lambda_bootstrap(j)).^3); %the PSD number function from Kok, PNAS, 2011

    %performing the lifetime resampling
    found_max_likelihood_T = false;
    while (~found_max_likelihood_T)
        try
            sample_array_T(j,1:no_T_data_sets) = randi([1 no_T_data_sets],1,no_T_data_sets);
            intercept = T_intercept(sample_array_T(j,:)); slope = T_slope(sample_array_T(j,:)); SE_intercept = T_SE_intercept(sample_array_T(j,:)); SE_slope = T_SE_slope(sample_array_T(j,:)); rho = T_rho(sample_array_T(j,:)); %setting the global fitting parameters equal to those for this particular bootstrap for the lifetime, in order to do the maximum likelihood estimation
            X0 =[mean(intercept),mean(slope)]; %initial parameters vector
            
            %defining the variables
            var_a = intercept; SE_var_a = SE_intercept;
            var_b = slope; SE_var_b = SE_slope;
            
            %this is the likelihood function, from combining the functions written out in fit_PSD_and_lifetime_functions; taking 1/L because fminunc finds the minimum, and we need the maximum
            L = @(X)1/prod((1./(2*pi*sqrt((SE_var_a.^2+tau_intercept^2).*(SE_var_b.^2+tau_slope^2).*(1-((rho.*SE_var_a.*SE_var_b+T_syst_error_corr*tau_intercept*tau_slope)./(sqrt(SE_var_a.^2+tau_intercept^2).*sqrt(SE_var_b.^2+tau_slope^2))).^2)))).*exp(-((var_a-X(1)).^2./(SE_var_a.^2+tau_intercept^2)+(var_b-X(2)).^2./(SE_var_b.^2+tau_slope^2)-2.*((rho.*SE_var_a.*SE_var_b+T_syst_error_corr*tau_intercept*tau_slope)./(sqrt(SE_var_a.^2+tau_intercept^2).*sqrt(SE_var_b.^2+tau_slope^2))).*(var_a-X(1)).*(var_b-X(2))./(sqrt(SE_var_a.^2+tau_intercept^2).*sqrt(SE_var_b.^2+tau_slope^2)))./(2.*(1-((rho.*SE_var_a.*SE_var_b+T_syst_error_corr*tau_intercept*tau_slope)./(sqrt(SE_var_a.^2+tau_intercept^2).*sqrt(SE_var_b.^2+tau_slope^2))).^2))));
            
            [X_min_T(j,1:2),fval(j),exitflag(j)] = fminunc(L,X0,options); %finding the parameters that maximize the lifetime likelihood function
            T0(j) = exp(X_min_T(j,1)); D_dep(j) = -1/X_min_T(j,2); %obtaining the fitting parameters
        if (fval(j)>1e-25)
            found_max_likelihood_T = true;
        else
            no_failed_max_likelihood_T = no_failed_max_likelihood_T+1;
            'problem in lifetime maximum likelihood method: fval too low'
            failed_sample_array_T(no_failed_max_likelihood_T,1:no_T_data_sets) = sample_array_T(j,1:no_T_data_sets);
        end
        catch
            no_failed_max_likelihood_T = no_failed_max_likelihood_T+1;
            'error in lifetime maximum likelihood method'
            failed_sample_array_T(no_failed_max_likelihood_T,1:no_T_data_sets) = sample_array_T(j,1:no_T_data_sets);
        end
    end
    T(j,:) = (T0(j)/daysInYears).*exp(-D./D_dep(j)); %the dust lifetime in years
    
    %determining the integrals in L_tot and F_tot
    prod_func_denominator = PSD_V(j,:).*(Qe(j,:)./D).*T(j,:);
    denominator(j) = 10^6*trapz(D,prod_func_denominator); %integrating to get the denominator, with units of years/m
    denominator_D10(j) = 10^6*trapz(D(1:i_D10),prod_func_denominator(1:i_D10)); %integrating to get the denominator, with units of years/m; PM10 only
    F_eff(j)=(10^(-9))*(2*rho_d(j)*A_Earth/3)/denominator(j); %global dust emission rate per unit of optical depth, in Tg/year
    
    prod_func_Ltot = PSD_V(j,:).*T(j,:);
    integral_Ltot(j) = trapz(D,prod_func_Ltot);
    integral_Ltot_D10(j) = trapz(D(1:i_D10),prod_func_Ltot(1:i_D10)); %PM10 only
    L_eff(j)=(10^(-9))*(2*rho_d(j)*A_Earth/3)*integral_Ltot(j)/denominator(j); %global dust loading per unit of optical depth, in Tg
    
    %calculating the size-resolved emission flux and atmospheric load
    tau_d_bootstrap(j) = DAOD(ceil(size(DAOD,1)*rand(1))); %choosing the dust AOD from the pdf from Ridley et al. (2016)
    F_emit(j) = F_eff(j)*tau_d_bootstrap(j); %the total dust flux from the bootstrap, in Tg/year
    L_atm(j) = tau_d_bootstrap(j)*L_eff(j); %the estimated atmospheric dust loading in Tg
    column_burden(j) = 1e12*L_atm(j)/A_Earth; %the estimated column burden in g/m2
    MEE(j) = tau_d_bootstrap(j)/column_burden(j); %the dust mass extinction efficiency
    PSD_lnV_emit(j,1:size(D,2)) = D.*F_emit(j).*PSD_V(j,:); %the bootstrap dust size distribution (not normalized, dV_emit/dlnD) at emission
    PSD_lnV_emit_norm(j,1:size(D,2)) = D.*PSD_V(j,:); %the *normalized* bootstrap dust size distribution (dV_emit/dlnD) at emission
    PSD_V_emit_fct = @(D) F_emit(j).*(1./Cv_bootstrap(j)).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %emitted size distribution function, to be used to calculate fraction in each bar chart bin
    PSD_lnV_load(j,1:size(D,2)) = PSD_lnV_emit(j,1:size(D,2)).*T(j,:); %the bootstrap dust volume size distribution (not normalized, dV_load/dlnD) in the atmosphere
    PSD_lnV_load_norm(j,1:size(D,2)) = (F_emit(j)/L_atm(j))*PSD_lnV_emit_norm(j,1:size(D,2)).*T(j,:); %the normalized bootstrap dust size distribution in the atmosphere
    PSD_lnN_load(j,1:size(D,2)) = 1e9*PSD_lnV_load(j,1:size(D,2))./((pi/6)*rho_d(j).*(1e-6*D).^3); %the bootstrap dust number size distribution (not normalized, dN_load/dlnD) in the atmosphere
    PSD_V_load_fct = @(D) (T0(j)/daysInYears).*exp(-D./D_dep(j)).*F_emit(j).*(1./Cv_bootstrap(j)).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %atmospheric load size distribution function, to be used to calculate fraction in each bar chart bin, in Tg
    PSD_N_load_fct = @(D) (1e9./((pi/6)*rho_d(j).*(1e-6*D).^3)).*(T0(j)/daysInYears).*exp(-D./D_dep(j)).*F_emit(j).*(1./Cv_bootstrap(j)).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %atmospheric load size distribution function, to be used to calculate fraction in each bar chart bin; 1e9 is there to convert Tg to kg
    PSD_lnV_AOD(j,1:size(D,2)) = 1e9*PSD_lnV_load(j,1:size(D,2)).*Qe(j,:)*3./(2*rho_d(j)*D*1e-6*A_Earth); %the bootstrap size distribution contribution to AOD
    PSD_lnV_AOD_norm(j,1:size(D,2)) = (1/L_atm(j))*1e9*PSD_lnV_load(j,1:size(D,2)).*Qe(j,:)*3./(2*rho_d(j)*D*1e-6*A_Earth); %the bootstrap size distribution contribution to AOD, per Tg loading
    PSD_V_AOD_fct = @(D) (1e9.*Qe(j,D)*3./(2*rho_d(j)*D*1e-6*A_Earth))*(T0(j)/daysInYears).*exp(-D./D_dep(j)).*F_emit(j).*(1./Cv_bootstrap(j)).*(1+erf(log(D./Ds_bootstrap(j))/(sqrt(2)*log(sigma_s_bootstrap(j))))).*exp(-(D/lambda_bootstrap(j)).^3); %atmospheric load size distribution function, to be used to calculate fraction in each bar chart bin
    %calculating the relative and cumulative contributions of dust of different sizes to global DAOD
    AOD_rel_contribution(j,:) = PSD_V_emit_fct(D).*(Qe(j,:)./D).*T(j,:)/trapz(D,PSD_V_emit_fct(D).*(Qe(j,:)./D).*T(j,:));
    AOD_rel_contribution_cum(j,1) = 0;
    for p=2:size(D,2)
        AOD_rel_contribution_cum(j,p) = trapz(D(1:p),AOD_rel_contribution(j,1:p));
    end    
    scattering_enhancement (j) = sum((PSD_lnV_load(j,:)./D).*(1./D).*Qe(j,:))/sum((PSD_lnV_load(j,:)./D).*(1./D).*Qe_Mie); %the enhancement in scattering over the spherical case, for particles with the same volume
    %fraction of emission rate, load, and AOD due to clay:
    for p=1:size(D_min,2) %calculating the content of each bin for the bar charts
        bin_emit(j,p) = integral(PSD_V_emit_fct,D_min(p),D_max(p));
        bin_load(j,p) = integral(PSD_V_load_fct,D_min(p),D_max(p));
        N_bin_load(j,p) = integral(PSD_N_load_fct,D_min(p),D_max(p));
        m_part(j,p) = (pi/6)*rho_d(j)*(1e-12*D_min(p)*D_max(p))^(3/2);
        N_bin_load_alt(j,p) = 1e9*bin_load(j,p)/m_part(j,p);
        bin_AOD(j,p) = trapz(D(D_i_bin_start(p):D_i_bin_end(p)),PSD_lnV_AOD(j,D_i_bin_start(p):D_i_bin_end(p))./D(D_i_bin_start(p):D_i_bin_end(p))); %cannot use 'integral' because Qe is not a function. Only trust this when accurate is turned on.
    end
    clay_emit_fr(j) = sum(bin_emit(j,1:4))/sum(bin_emit(j,:));
    clay_load_fr(j) = sum(bin_load(j,1:4))/sum(bin_load(j,:));
    clay_N_load_fr(j) = sum(N_bin_load(j,1:4))/sum(N_bin_load(j,:));    
    clay_AOD_fr(j) = sum(bin_AOD(j,1:4))/sum(bin_AOD(j,:));
    F_emit_D10(j) = sum(bin_emit(j,1:find(D_max==10))); %Emission of PM10 dust only
    L_atm_D10(j) = sum(bin_load(j,1:find(D_max==10))); %load of PM10 dust only
    column_burden_D10(j) = 1e12*L_atm_D10(j)/A_Earth; %the estimated column burden in g/m2
    MEE_D10(j) = tau_d_bootstrap(j)/column_burden_D10(j); %the dust mass extinction efficiency
    
    L_atm_PM1_geom(j) = integral(PSD_V_load_fct,min(D),1);
    L_atm_PM25_geom(j) = integral(PSD_V_load_fct,min(D),2.5);
    for k=1:no_DRE_models %calculating the AOD in each bin of each DRE model
        for p=1:no_DRE_bins(k)
            if (D_DRE_upper(k,p)<20)
                DRE_bin_load(k,j,p) = integral(PSD_V_load_fct,D_DRE_lower(k,p),D_DRE_upper(k,p));
                DRE_bin_AOD(k,j,p) = trapz(D(D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p)),PSD_lnV_AOD(j,D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p))./D(D_DRE_i_bin_start(k,p):D_DRE_i_bin_end(k,p))); %cannot use 'integral' because Qe is not a function. Only trust this when accurate is turned on.
            end
            if (p==no_DRE_bins(k) && D_DRE_upper(k,p)<20) %the last bin does not extend to 20 um
                DRE_bin_AOD_PM20(k,j) = trapz(D(D_DRE_i_bin_end(k,no_DRE_bins(k)):end),PSD_lnV_AOD(j,D_DRE_i_bin_end(k,no_DRE_bins(k)):end)./D(D_DRE_i_bin_end(k,no_DRE_bins(k)):end)); %DAOD due to the PM20 size range not covered in CESM
            else %the last bin does extend to 20 um
                DRE_bin_AOD_PM20(k,j) = trapz(D(D_DRE_i_bin_end(k,no_DRE_bins(k)-1):end),PSD_lnV_AOD(j,D_DRE_i_bin_end(k,no_DRE_bins(k)-1):end)./D(D_DRE_i_bin_end(k,no_DRE_bins(k)-1):end)); %DAOD due to the PM20 size range not covered in GISS
            end
        end
    end %for, calculating the AOD in each bin of each DRE model
end %for, looping over the resampling iterations

%saving the AOD for each model particle bin:
save('Model_AOD_per_bin.mat','DRE_bin_AOD','DRE_bin_AOD_PM20','rho_d','tau_d_bootstrap');
save('Vars_for_AeroCom_to_DRE_model_AOD.mat','D','D_m','Qe','Qe_Mie','rho_d'); %the variables needed to calculate the DAOD from AeroCom models, mapped onto the DRE models
if (no_failed_max_likelihood_T>0)
    warning('max likelihood method had failures for lifetime')
    no_failed_max_likelihood_T
end

%calculating the atmospheric dust loading and its error
A=sort(L_atm);
L_atm_avg = A(round(0.5*j_end));
L_atm_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
L_atm_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
L_atm_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
L_atm_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma
%PM10 (geometric) only
A=sort(L_atm_D10); 
L_atm_D10_avg = A(round(0.5*j_end));
L_atm_D10_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
L_atm_D10_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
L_atm_D10_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
L_atm_D10_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma

%PM2.5 (geometric) only
A=sort(L_atm_PM25_geom); 
L_atm_PM25_geom_avg = A(round(0.5*j_end));
L_atm_PM25_geom_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
L_atm_PM25_geom_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
L_atm_PM25_geom_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
L_atm_PM25_geom_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma

%calculating the global dust loading per unit of optical depth, in Tg, and its error
A=sort(L_eff);
L_eff_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
L_eff_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
L_eff_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
L_eff_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
L_eff_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma

%calculating the optical depth per Tg of dust loading and its error
DAOD_per_Tg = 1./L_eff;
A=sort(DAOD_per_Tg);
DAOD_per_Tg_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
DAOD_per_Tg_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
DAOD_per_Tg_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
DAOD_per_Tg_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
DAOD_per_Tg_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma

%calculating the mass extinction efficiency and its error
A=sort(MEE);
MEE_avg = A(round(0.5*j_end));
MEE_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
MEE_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
MEE_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
MEE_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma
%PM10 (geometric) only
A=sort(MEE_D10); 
MEE_D10_avg = A(round(0.5*j_end));
MEE_D10_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
MEE_D10_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
MEE_D10_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
MEE_D10_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma

%calculating the global emission rate and its error
A=sort(F_emit);
F_emit_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
F_emit_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
F_emit_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
F_emit_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
F_emit_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma
%PM10 only
A=sort(F_emit_D10);
F_emit_D10_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
F_emit_D10_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
F_emit_D10_neg_1sigma = A(max(1,(floor(0.5*(1+erf(-1/sqrt(2)))*j_end)))); %lower uncertainty value at -1 sigma
F_emit_D10_pos_1sigma = A(max(1,(ceil(0.5*(1+erf(1/sqrt(2)))*j_end)))); %upper uncertainty value at +1 sigma
F_emit_D10_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma

F_eff_rel_err = std(F_eff)/mean(F_eff); %relative error in F_eff
F_tot = mean(DAOD)*mean(F_eff); %best estimate of total dust emission rate
F_std = F_tot*sqrt((std(DAOD)/mean(DAOD))^2+F_eff_rel_err^2); %error in estimate of total dust emission rate

L_eff_rel_err = std(L_eff)/mean(L_eff); %relative error in L_eff
L_tot = mean(DAOD)*mean(L_eff); %best estimate of total dust emission rate
L_std = L_tot*sqrt((std(DAOD)/mean(DAOD))^2+L_eff_rel_err^2); %error in estimate of total dust emission rate

T_mean = (mean(T0(j))/daysInYears).*exp(-D./mean(D_dep)); %the dust lifetime in years

A=sort(clay_emit_fr);
clay_emit_fr_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
clay_emit_fr_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
clay_emit_fr_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma

A=sort(clay_load_fr);
clay_load_fr_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
clay_load_fr_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
clay_load_fr_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma

A=sort(clay_N_load_fr);
clay_N_load_fr_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
clay_N_load_fr_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
clay_N_load_fr_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma

A=sort(clay_AOD_fr);
clay_AOD_fr_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
clay_AOD_fr_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
clay_AOD_fr_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma

%calculating extinction properties
Qe_mean=mean(Qe_TAMU);
Qe_err=std(Qe_TAMU);

A=sort(geometric_ext_efficiency);
geometric_ext_efficiency_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
geometric_ext_efficiency_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
geometric_ext_efficiency_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %upper uncertainty value at +2 sigma    

A=sort(scattering_enhancement);
scattering_enhancement_avg = A(round(0.5*j_end)); %lower uncertainty value at -1 sigma
scattering_enhancement_pos_2sigma = A(max(1,(ceil(0.5*(1+erf(2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma
scattering_enhancement_neg_2sigma = A(max(1,(floor(0.5*(1+erf(-2/sqrt(2)))*j_end)))); %lower uncertainty value at -2 sigma

write_data; %writes out data to files