function [T_intercept, T_slope, T_SE_intercept, T_SE_slope, T_rho, T0_max_likelihood, D_dep_max_likelihood, Cn, Ds, sigma_s, lambda, lambda_defined, PSD_covariance_matrix, Ds_max_likelihood, sigma_s_max_likelihood, lambda_max_likelihood] = fit_PSD_and_lifetime_functions (~)

%this function fits the lifetime and emitted dust PSD data sets using the
%statistical model described in the supplement to Kok et al. (2017),
%obtaining the MLE

global PM10_only eliminate_Yue_2010 calculate_likelihood_function_gradient

PM10_only = false; %sets whether to restrict the data to PM10 only; should run for PM20 because consistent with Kok (2011), and will produce stronger constraints on PM10
options = optimoptions(@fminunc,'GradObj','off','MaxFunEvals',2000); %sets whether to use gradient method for obtaining function minimum
differentialWeight = true; %sets whether to use differential weights in correcting PSD data to best fit function
eliminate_Yue_2010 = false; %if set to true, eliminates Yue et al. (2010) lifetime data, because (i) the model is not of the same quality as the other models, which all have extensive working groups and quality control measures, (ii) the study is older, and (iii) the study was run at very coarse resolution
calculate_likelihood_function_gradient = false; %sets whether to use the local gradient for the algorithm in the max likelihood method
D_end = 20; %upper limit on integration over particle size, in um
plot_data = true; %switch that sets whether to plot the data and the maximum likelihood estimate functions
if (plot_data)
    figure(1); clf;
    figure(2); clf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%      reading in the PSD data       %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PSD_x, PSD_y, PSD_y_err, N_PSD, no_PSD_data_sets] = load_PSD_data;
[T_x, T_y, T_x_st, T_x_en, N_T, no_T_data_sets] = load_lifetime_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtaining the most likely estimates for the lifetime function, from the maximum likelihood method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:no_T_data_sets
    X = T_x(i,1:N_T(i));
    Y = log(T_y(i,1:N_T(i)));
    [T_intercept(i), T_slope(i), T_SE_intercept(i), T_SE_slope(i), T_covariance(i)] = linearfit (X, Y, 0);    
    clear X; clear Y;    
end
T_rho = T_covariance.^2./(T_SE_intercept.*T_SE_slope);
intercept = T_intercept; slope = T_slope; SE_intercept = T_SE_intercept; SE_slope = T_SE_slope; rho = T_rho; %setting the global fitting parameters equal to those for the lifetime in order to do the maximum likelihood estimation
X0 =[mean(intercept),mean(slope),std(intercept),std(slope),mean(rho)]; %initial parameters vector

%defining the variables
var_a = intercept; SE_var_a = SE_intercept;
var_b = slope; SE_var_b = SE_slope;

%this is the likelihood function, from combining the functions above; taking 1/L because fminunc finds the minimum, and we need the maximum
L = @(X)1/prod((1./(2*pi*sqrt((SE_var_a.^2+X(3)^2).*(SE_var_b.^2+X(4)^2).*(1-((rho.*SE_var_a.*SE_var_b+X(5)*X(3)*X(4))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2))).^2)))).*exp(-((var_a-X(1)).^2./(SE_var_a.^2+X(3)^2)+(var_b-X(2)).^2./(SE_var_b.^2+X(4)^2)-2.*((rho.*SE_var_a.*SE_var_b+X(5)*X(3)*X(4))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2))).*(var_a-X(1)).*(var_b-X(2))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2)))./(2.*(1-((rho.*SE_var_a.*SE_var_b+X(5)*X(3)*X(4))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2))).^2))));

%finding the minimum parameters 
[X_min,fval] = fminunc(L,X0,options); %finding the parameters that maximize the lifetime likelihood function
intercept_max_likelihood = X_min(1); slope_max_likelihood = X_min (2);
tau_intercept = X_min(3); tau_slope = X_min (4);
T_syst_error_corr = X_min(5); %the correlation between systematic errors in intercept and slope
T0_max_likelihood = exp(X_min(1)); D_dep_max_likelihood = -1/X_min(2);
tau_T0 = T0_max_likelihood*X_min(3); tau_D_dep = (1/X_min(2)^2)*X_min(4); %from error propagation

figure(2); clf; semilogy(1,1); hold on; 
if (plot_data)
    for i=1:no_T_data_sets
        x_plot = 0:30;
        y_plot = exp(intercept(i))*exp(x_plot*slope(i));
        y_plot2 = T0_max_likelihood*exp(-x_plot/D_dep_max_likelihood);
        figure (1); clf; semilogy(T_x(i,1:N_T(i)),T_y(i,1:N_T(i)),'sb'); hold on; semilogy(x_plot,y_plot); 
    	figure (2); semilogy(T_x(i,1:N_T(i)),T_y(i,1:N_T(i)),'sb'); plot(x_plot,y_plot2);
    end
end
figure(2); plot(x_plot,y_plot2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating the maximum likelihood estimates for the dust PSD parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:5
    switch i
    case 1 %Gillette et al. (1972, 1974); Gillette (1974)
        Ds(i) = 1.86531;
        sigma_s(i) = 1.94868;
        lambda(i) = 100;
        Cn(i) = 1.84224;
        lambda_defined(i) = false;
        PSD_covariance_matrix(i,1:2,1:2) = [0.13561, 0.10748; 0.10748, 0.17412]; 
    case 2 %Fratini et al. (2007)
        Ds(i) = 1.62807;  %from Origin fit
        sigma_s(i) = 2.57463;  %from Origin fit
        lambda(i) = 100;  %from Origin fit
        Cn(i) = 2.00008;  %from Origin fit
        lambda_defined(i) = false;
        PSD_covariance_matrix(i,1:2,1:2) = [0.50566, 0.53978; 0.53978, 0.681];  %from Origin fit        
    case 3 %Sow et al. (2009)
        Ds(i) = 3.14381; 
        sigma_s(i) = 2.99276;
        lambda(i) = 12.95786;
        Cn(i) = 1.07453;
        lambda_defined(i) = true;
        PSD_covariance_matrix(i,1:3,1:3) = [0.36389, 0.21035, 0.08295; 0.21035, 0.14545, 0.05191; 0.08295, 0.05191, 0.40449];  %from Origin fit        
    case 4 %Shao et al. (2011)
        Ds(i) = 1.42472;  %from Origin fit
        sigma_s(i) = 1.81262;  %from Origin fit
        lambda(i) = 100;  %from Origin fit
        Cn(i) = 1.75509;  %from Origin fit
        lambda_defined(i) = false;
        PSD_covariance_matrix(i,1:2,1:2) = [0.13178, 0.14844; 0.14844, 0.24334];  %from Origin fit        
    case 5 %Rosenberg et al. (2014)
        Ds(i) = 0.60937; 
        sigma_s(i) = 2.16819;
        lambda(i) = 29.78776;
        Cn(i) = 1.88088;
        lambda_defined(i) = true;
        PSD_covariance_matrix(i,1:3,1:3) = [0.02614, -0.09235, -2.39769e-4; -0.09235, 1.41467, 0.00671; -2.39769e-4, 0.00671, 14.89117];        
    end %switch, i
    y_int = @(D) (1+erf(log(D./Ds(i))/(sqrt(2)*log(sigma_s(i))))).*exp(-(D./lambda(i)).^3); %the function to be integrated
end %for i
Ds_err = sqrt(squeeze(PSD_covariance_matrix(:,1,1))');
sigma_s_err = sqrt(squeeze(PSD_covariance_matrix(:,2,2))');
lambda_err = sqrt(squeeze(PSD_covariance_matrix(:,3,3))');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   calculating the best fit for the PSD   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first calculating lambda, by assuming it's drawn from a normal distribution, with systematic error tau_lambda added
lambda_fct = lambda(lambda_defined); SE_lambda_fct = lambda_err(lambda_defined);
X0 = [mean(lambda_fct), std(lambda_fct)];
f = @(X)1/(prod((1./sqrt(2*pi*(SE_lambda_fct.^2+X(2)^2))).*exp(-((lambda_fct-X(1)).^2)./(2*(SE_lambda_fct.^2+X(2)^2)))));
[lambda_min,lambda_fval,lambda_exitflag] = fminunc(f,X0, options); %finding the parameters that maximize the PSD likelihood function
lambda_min

lambda_max_likelihood = lambda_min(1); tau_lambda = lambda_min(2); %most likely value of lambda
%the below corrects lambda for those data sets for which it's undefined, and also corrects Cv accordingly
lambda(find(lambda_defined==false)) = lambda_max_likelihood; %setting the lambda for ambiguous data sets equal to the MLE estimate
lambda_err(find(lambda_defined==false)) = tau_lambda; %setting the error of lambda for ambiguous data sets equal to the systematic error

for i=1:5 %calculating Cv(i) for each study
	y_int = @(D) (1+erf(log(D./Ds(i))/(sqrt(2)*log(sigma_s(i))))).*exp(-(D./lambda(i)).^3); 
    Cv(i) = integral(y_int,0.2,D_end); 
end

%then calculating sigma_s and Ds, which are correlated and thus described by a bivariate normal distribution
PSD_rho = squeeze(PSD_covariance_matrix(:,1,2))'./(Ds_err.*sigma_s_err);
rho = PSD_rho;
A = corrcoef(Ds,sigma_s);
rho_syst_estimate = A(1,2); %taking the correlation between the sigma_s and Ds systematic errors as the correlation between sigma_s and Ds
X0 =[mean(Ds),mean(sigma_s),std(Ds),std(sigma_s)]; %initial parameters vector

%defining the variables
var_a = Ds; SE_var_a = Ds_err;
var_b = sigma_s; SE_var_b = sigma_s_err;

%this is the likelihood function, from combining the functions above; taking 1/L because fminunc finds the minimum, and we need the maximum
L = @(X)1/prod((1./(2*pi*sqrt((SE_var_a.^2+X(3)^2).*(SE_var_b.^2+X(4)^2).*(1-((rho.*SE_var_a.*SE_var_b+rho_syst_estimate*X(3)*X(4))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2))).^2)))).*exp(-((var_a-X(1)).^2./(SE_var_a.^2+X(3)^2)+(var_b-X(2)).^2./(SE_var_b.^2+X(4)^2)-2.*((rho.*SE_var_a.*SE_var_b+rho_syst_estimate*X(3)*X(4))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2))).*(var_a-X(1)).*(var_b-X(2))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2)))./(2.*(1-((rho.*SE_var_a.*SE_var_b+rho_syst_estimate*X(3)*X(4))./(sqrt(SE_var_a.^2+X(3)^2).*sqrt(SE_var_b.^2+X(4)^2))).^2))));

%finding the minimum parameters from minimizing the likelihood function
[soilPSD_X_min,soilPSD_fval,soilPSD_exitflag]  = fminunc(L,X0,options); %finding the parameters that maximize the lifetime likelihood function
Ds_max_likelihood = soilPSD_X_min(1); sigma_s_max_likelihood = soilPSD_X_min(2); tau_Ds = soilPSD_X_min(3); tau_sigma_s = soilPSD_X_min(4);

%calculating Cn_max_likelihood and Cv_max_likelihood
y_int = @(D) (1./D.^3).*(1+erf(log(D./Ds_max_likelihood)/(sqrt(2)*log(sigma_s_max_likelihood)))).*exp(-(D/lambda_max_likelihood).^3); %the function to be integrated
Cn_max_likelihood = integral(y_int,0.2,D_end); %calculating the normalization constant for the number PSD
y_int = @(D) (1+erf(log(D./Ds_max_likelihood)/(sqrt(2)*log(sigma_s_max_likelihood)))).*exp(-(D/lambda_max_likelihood).^3); %the function to be integrated
Cv_max_likelihood = integral(y_int,0.2,D_end); %calculating the normalization constant for the volume PSD


%obtain the correction factor by which the PSD data needs to be multiplied to maximize agreement against the most likely normalized PSD
PSD_weight = null(1);
D = [logspace(-1, log10(10), 101),logspace(1, log10(20), 16)]; %particle diameter in um
for i = 1:no_PSD_data_sets
    PSD_weight(1+size(PSD_weight,2):N_PSD(i)+size(PSD_weight,2)) = (1/N_PSD(i))*ones(1,N_PSD(i));
    PSD_lnV_max_likelihood = (D./Cv(i)).*(1+erf(log(D./Ds(i))/(sqrt(2)*log(sigma_s(i))))).*exp(-(D./lambda(i)).^3); %the the maximum likelihood PSD volume function for this particular data set, from Kok, PNAS, 2011
    [PSD_corr_fact(i),PSD_y_corr(i,1:N_PSD(i)), PSD_y_err_corr(i,1:N_PSD(i))] = calc_PSD_corr_fact(PSD_x(i,1:N_PSD(i)), PSD_y(i,1:N_PSD(i)), PSD_y_err(i,1:N_PSD(i)), PSD_weight, D, PSD_lnV_max_likelihood, 1);
    figure(1); clf; loglog(PSD_x(i,1:N_PSD(i)),PSD_y_corr(i,1:N_PSD(i)),'ok'); hold on; loglog(D,PSD_lnV_max_likelihood,'-k');
end
axis([0.2 20 2e-4 2]);

%writing out normalized PSD data; need to split up first (Gillette) data set into its three constituent data sets
fid = fopen('PSD_data_normalized.txt','wt');
PSD_x_write(4:7,1:size(PSD_x,2))= PSD_x(2:5,:); PSD_y_write(4:7,1:size(PSD_x,2))= PSD_y_corr(2:5,:); PSD_y_err_write(4:7,1:size(PSD_x,2))= PSD_y_err_corr(2:5,:);
N_PSD_write(4:7) = N_PSD(2:5);
%Gilette et al. (1972)
PSD_x_write(1,1:4) = PSD_x(1,1:4); PSD_y_write(1,1:4) = PSD_y_corr(1,1:4); PSD_y_err_write(1,1:4) = PSD_y_err_corr(1,1:4);
N_PSD_write(1) = 4;
%Gillette et al. (1974)
PSD_x_write(2,1:4) = PSD_x(1,5:8); PSD_y_write(2,1:4) = PSD_y_corr(1,5:8); PSD_y_err_write(2,1:4) = PSD_y_err_corr(1,5:8);
N_PSD_write(2) = 4;
%Gillette (1974)
N_PSD_write(3) = 4;
PSD_x_write(3,1:4) = PSD_x(1,9:12); PSD_y_write(3,1:4) = PSD_y_corr(1,9:12); PSD_y_err_write(3,1:4) = PSD_y_err_corr(1,9:12);

m_array = [1, 2, 3, 5, 6, 4, 7]; %setting the order of data sets to be written out
N_max = max(N_PSD_write);
for p=1:N_max
    for i=1:size(m_array,2)
        m=m_array(i);
        N=N_PSD_write(m);
        if (p<=N_PSD_write(m))
            fprintf(fid,'%1.4e %1.4e %1.4e ', PSD_x_write(m,p), PSD_y_write(m,p), PSD_y_err_write(m,p));
        else
            fprintf(fid,'. . . '); %no data
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);

%saving the data to a datafile, to be read in by the main module
filename = 'lifetime_and_PSD_fits.mat';
save(filename,'T_intercept', 'T_slope', 'T_SE_intercept', 'T_SE_slope', 'T_rho', 'T0_max_likelihood', 'D_dep_max_likelihood', 'tau_intercept', 'tau_slope','T_syst_error_corr','Cn', 'Ds', 'Ds_err', 'sigma_s', 'sigma_s_err', 'lambda', 'lambda_err', 'lambda_defined', 'PSD_covariance_matrix', 'Ds_max_likelihood', 'sigma_s_max_likelihood', 'lambda_max_likelihood', 'Cv_max_likelihood','rho_syst_estimate','tau_Ds','tau_sigma_s');

%plotting the data
x_plot = min(min(PSD_x(PSD_x>0))):0.2:max(max(PSD_x(PSD_x>0)));
figure(1); clf; loglog(1,1); hold on;
figure(2); clf; loglog(1,1); hold on;
i_start = 1;
for i=1:no_PSD_data_sets
    i_end = i_start + N_PSD(i) - 1;
    PSD_min = @(D) (D./Cv(i)).*(1+erf(log(D./Ds(i))/(sqrt(2)*log(sigma_s(i))))).*exp(-(D/lambda(i)).^3); %the PSD function from Kok, PNAS, 2011
    figure(2); loglog(x_plot, PSD_min(x_plot),'k'); hold on;
    X = PSD_x(i,1:N_PSD(i));
    Y = PSD_y_corr(i,1:N_PSD(i));
    %Y_err = PSD_y_err(i,1:N_PSD(i));
    figure(1); loglog(X,Y,'rs');    
    figure(2); loglog(X,Y,'rs');
    i_start = i_end + 1;
    PSD_min = @(D) (D./Cv(i)).*(1+erf(log(D./Ds(i))/(sqrt(2)*log(sigma_s(i))))).*exp(-(D/lambda(i)).^3); %the PSD function from Kok, PNAS, 2011
    figure(2); loglog(x_plot, PSD_min(x_plot),'k'); hold on;        
end
PSD_MSE = @(D) (D./Cv_max_likelihood).*(1+erf(log(D./Ds_max_likelihood)/(sqrt(2)*log(sigma_s_max_likelihood)))).*exp(-(D/lambda_max_likelihood).^3); %the PSD function from Kok, PNAS, 2011figure(1); loglog(x_plot, PSD_min(x_plot),'k'); hold on;
figure(1); loglog(x_plot, PSD_MSE(x_plot),'k');
axis([0.9*min(min(PSD_x(PSD_x>0))) 1.1*max(max(PSD_x(PSD_x>0))) 0.9*min(min(PSD_y(PSD_y>0))) 1.1*max(max(PSD_y(PSD_y>0)))]);