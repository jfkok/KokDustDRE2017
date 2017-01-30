function [x, y, x_st, x_en, N, no_data_sets] = load_lifetime_data

%function that loads the model data on the dust lifetime

global PM10_only

%loading the lifetime data
%1 = ModelE, 2 = CESM, 3 = GEOS-Chem, 4 = WRF-Chem, 5 = HadGEM, 6 = MERRAero, 7 = MOZART, 8 = UMI, 9 = GMOD

%from Perlwitz et al. (2014, ACPD, part 1), Fig. S12, AMF method
m = 1; %first data set
if (PM10_only)
    N(m) = 4; %number of data points
    x(m,1:N(m)) = [0.63246, 2.82843, 5.65685, 11.31371]; %Dmean
    y(m,1:N(m)) = [5.83, 5.18, 3.7, 1.42]; %lifetime in days
    x_st(m,1:N(m)) = [0.1, 2, 4, 8]; %bin lower limit
    x_en(m,1:N(m)) = [2, 4, 8, 16]; %bin upper limit    
    L = [0.71, 2.12, 2.74, 2.44]; %from Fig. 7 in Perlwitz et al. (2014, ACPD, part 1)
    emitted_flux = L./(y(m,1:N(m))/365); %flux(kg/year) = load (kg) / lifetime (years)
else
    N(m) = 5; %number of data points    
    x(m,1:N(m)) = [0.63246, 2.82843, 5.65685, 11.31371, 22.62742]; %Dmean
    y(m,1:N(m)) = [5.83, 5.18, 3.7, 1.42, 0.26]; %lifetime in days
    x_st(m,1:N(m)) = [0.1, 2, 4, 8, 16]; %bin lower limit
    x_en(m,1:N(m)) = [2, 4, 8, 16, 32]; %bin upper limit
    L = [0.71, 2.12, 2.74, 2.44, 0.81]; %from Fig. 7 in Perlwitz et al. (2014, ACPD, part 1)
    emitted_flux = L./(y(m,1:N(m))/365); %flux(kg/year) = load (kg) / lifetime (years)
end %if, PM10_only

%from Kok et al. (2014)
m = 2; %mth data set
N(m) = 4;
x(m,1:N(m)) = [0.31623, 1.58114, 3.53553, 7.07107]; %Dmean
y(m,1:N(m)) = [9.95, 9.35, 5.46, 1.94]; %lifetime in days
x_st(m,1:N(m)) = [0.1, 1, 2.5, 5]; %bin lower limit
x_en(m,1:N(m)) = [1, 2.5, 5, 10]; %bin upper limit
emitted_flux = 1e12*[0.0324    0.2563    0.8162    1.8415]; %emitted flux in kg

%from Ridley (personal communication, 2014)
m = 3; %mth data set
N(m) = 4;    
x(m,1:N(m)) = [0.63246, 2.68328, 4.64758, 8.48528]; %Dmean
y(m,1:N(m)) = [7.15603, 6.3739, 4.34741, 1.12857]; %lifetime in days
x_st(m,1:N(m)) = [0.2, 2, 3.6, 6]; %bin lower limit
x_en(m,1:N(m)) = [2, 3.6, 6, 12]; %bin upper limit
emitted_flux = [160.3358, 331.794, 422.896, 396.135];

%from Chun Zhao, (personal communication, 2015)
m = 4; %mth data set
N(m) = 4;    
x_st(m,1:N(m)) = [0.039, 0.625, 2.5, 5]; %lower bin limits in um, from Zhao et al., 2010; for deposition, bins 1-4 and 5-6 are combined
x_en(m,1:N(m)) = [0.625, 2.5, 5, 10]; %in um
x(m,1:N(m)) = sqrt(x_st(m,1:N(m)).*x_en(m,1:N(m))); %in um
y(m,1:N(m)) = [9.6224    6.2381    3.2246    2.3690]; %lifetime in days
emitted_flux = 1e9*[17.8  412.9  1217.7  1970.8]; %emitted flux in kg

%from Karsten's HadGEM run
m = 5; %mth data set    
N(m) = 5;
x(m,1:N(m)) = [0.11243, 0.35553, 1.12428, 3.55528, 11.24278]; %Dmean
y(m,1:N(m)) = [14.67, 14.977, 12.101, 6.91, 1.01767]; %lifetime in days 
x_st(m,1:N(m)) = [0.0632, 0.2, 0.632, 2, 6.32]; %bin lower limit
x_en(m,1:N(m)) = [0.2, 0.632, 2, 6.32, 20]; %bin upper limit    
    
%from MERRAero
m = 6; %mth data set
if (PM10_only)
    N(m) = 4;
    x_st(m,1:N(m)) = [0.2, 2, 3.6, 6]; %lower bin limits in um, from Zhao et al., 2010; for deposition, bins 1-4 and 5-6 are combined
    x_en(m,1:N(m)) = [2, 3.6, 6, 12]; %in um
    x(m,1:N(m)) = [1.272, 2.64, 4.6, 8.34]; %in um
    y(m,1:N(m)) = [10.926, 9.506, 6.592, 2.354]; %lifetime in days
else
    N(m) = 5;
    x_st(m,1:N(m)) = [0.2, 2, 3.6, 6, 12]; %lower bin limits in um, from Zhao et al., 2010; for deposition, bins 1-4 and 5-6 are combined
    x_en(m,1:N(m)) = [2, 3.6, 6, 12, 20]; %in um
    x(m,1:N(m)) = [1.272, 2.64, 4.6, 8.34, 15.34]; %in um
    y(m,1:N(m)) = [10.926, 9.506, 6.592, 2.354, 0.484]; %lifetime in days
end

%MOZART, from Liu et al., Atm Env, 2009, table 1
m = 7; %mth data set
if (PM10_only)
    N(m) = 4;
    x_st(m,1:N(m)) = [0.2, 2, 3.6, 6]; %lower bin limits in um, from Zhao et al., 2010; for deposition, bins 1-4 and 5-6 are combined
    x_en(m,1:N(m)) = [2, 3.6, 6, 12]; %in um
    x(m,1:N(m)) = [0.632, 2.68, 4.65, 8.49]; %in um
    y(m,1:N(m)) = [17, 11, 4.8, 1.3]; %lifetime in days
else
    N(m) = 5;
    x_st(m,1:N(m)) = [0.2, 2, 3.6, 6, 12]; %lower bin limits in um, from Zhao et al., 2010; for deposition, bins 1-4 and 5-6 are combined
    x_en(m,1:N(m)) = [2, 3.6, 6, 12, 20]; %in um
    x(m,1:N(m)) = [0.632, 2.68, 4.65, 8.49, 15.49]; %in um
    y(m,1:N(m)) = [17, 11, 4.8, 1.3, 0.4]; %lifetime in days
end

%UMI, from Liu and Penner (2007)
m = 8; %mth data set
N(m) = 4;
x_st(m,1:N(m)) = [0.05, 0.63, 1.25, 2.5]; %lower bin limits in um, from Zhao et al., 2010; for deposition, bins 1-4 and 5-6 are combined
x_en(m,1:N(m)) = [0.63, 1.25, 2.5, 10]; %in um
x(m,1:N(m)) = [0.177, 0.89, 1.77, 5]; %in um
y(m,1:N(m)) = [7.726, 7.53, 8.34, 2.09]; %lifetime in days

    %from Yue et al. (2010)
    m = 9; %second data set    
    if (PM10_only)
        N(m) = 3;
        x(m,1:N(m)) = [0.78256, 2.99826, 6.80827]; %Deff
        x_st(m,1:N(m)) = [0.2, 2, 5]; %bin lower limit
        x_en(m,1:N(m)) = [2, 5, 10]; %bin upper limit
        y(m,1:N(m)) = [20.9, 14.5, 6.6]; %lifetime in days
        emitted_flux = [97, 213, 600]; %emitted flux in Tg, from Table 2 in Yue et al. (2009)
    else
        N(m) = 4;
        x(m,1:N(m)) = [0.78256, 2.99826, 6.80827, 13.578]; %Deff
        x_st(m,1:N(m)) = [0.2, 2, 5, 10]; %bin lower limit
        x_en(m,1:N(m)) = [2, 5, 10, 20]; %bin upper limit
        y(m,1:N(m)) = [20.9, 14.5, 6.6, 1.1]; %lifetime in days, from Table 2 in Yue et al., 2009
        emitted_flux = [97, 213, 600, 1025]; %emitted flux in Tg, from Table 2 in Yue et al. (2009)
    end  

no_data_sets = m; %number of model lifetime results