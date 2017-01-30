function [x, y, y_err, N, no_data_sets] = load_PSD_data

%function that loads the experimental data on the emitted dust PSD

global PM10_only 

%fitting a number size distribution
    %Gillette et al. (1972), (1974a), (1974b)
    m = 1; %second data set
    if (PM10_only)
        N(m) = 10;
        x(m,1:N(m)) = [0.82782, 1.89436, 4.58479, 9.86682, 1.72856, 3.54304, 6.90892, 1.74187, 3.3827, 6.77538];
        y(m,1:N(m)) = [0.19211, 0.10913, 0.0482, 0.01321, 0.16845, 0.09188, 0.02483, 0.15301, 0.10005, 0.0248];
        y_err(m,1:N(m)) = [0.10658, 0.06055, 0.02674, 0.00733, 0.06871, 0.00712, 0.00219, 0.02676, 0.01102, 0.00264];
    else
        N(m)=12;
        x(m,1:N(m)) = [0.82782, 1.89436, 4.58479, 9.86682, 1.72856, 3.54304, 6.90892, 13.72303, 1.74187, 3.3827, 6.77538, 13.54969];
        y(m,1:N(m)) = [0.19211, 0.10913, 0.0482, 0.01321, 0.16845, 0.09188, 0.02483, 0.01093, 0.15301, 0.10005, 0.0248, 0.00152];
        y_err(m,1:N(m)) = [0.10658, 0.06055, 0.02674, 0.00733, 0.06871, 0.00712, 0.00219, 0.00433, 0.02676, 0.01102, 0.00264, 2.57002E-4];
    end
    
    %Fratini et al. (2007)
    m = 2; %fourth data set
    N(m)=8;    
    %original Fratini data
    Frat_D = [0.27982, 0.31946, 0.37931, 0.42881, 0.46841, 0.51925, 0.64753, 0.75022, 0.84813, 0.9479, 1.24766, 1.74464, 2.2517, 2.74896, 3.50181, 4.49747, 5.49965, 6.49808];
    Frat_dNdlnD = [2.46314, 0.20156, 0.23409, 0.18418, 0.31462, 0.30847, 0.56179, 0.59304, 0.31382, 0.60805, 0.0679, 0.16125, 0.58143, 0.17774, 0.04283, 0.01974, 0.00616, 0.00167];
    Frat_dNdlnD_SE = [0.07788, 0.04174, 0.07375, 0.02163, 0.06492, 0.12082, 0.01726, 0.01487, 0.01428, 0.06261, 0.0042, 0.01184, 0.03785, 0.00668, 0.00413, 0.00347, 0.00131, 5.89733E-4];
    for i=1:8 %calculating the average of two adjacent bins; stripping the last bin because affected by cut-off size of inlet (personal communication with Fratini)
        Frat_D_2bin(i) = sqrt(Frat_D(2*i-1)*Frat_D(2*i));
        Frat_dNdlnD_2bin(i) = sqrt(Frat_dNdlnD(2*i-1)*Frat_dNdlnD(2*i));
        Frat_dNdlnD_SE_2bin(i) = sqrt(Frat_dNdlnD_SE(2*i-1)*Frat_dNdlnD_SE(2*i));
    end    
    x(m,1:N(m)) = Frat_D_2bin;
    y(m,1:N(m)) = Frat_dNdlnD_2bin;
    y_err(m,1:N(m)) = Frat_dNdlnD_SE_2bin;
        
    %Sow et al. (2009)
    m = 3; %fifth data set
    if (PM10_only)        
        N(m)=12;
        x(m,1:N(m)) = [0.34503, 0.44414, 0.57172, 0.72714, 0.88671, 1.27186, 1.78096, 2.43457, 3.47112, 4.4414, 6.14486, 8.70859];
        y(m,1:N(m)) = [0.70435, 0.23513, 0.23913, 0.23837, 0.24489, 0.23233, 0.21563, 0.16968, 0.104, 0.0711, 0.03218, 0.01289];
        y_err(m,1:N(m)) = [0.28337, 0.03855, 0.05773, 0.05338, 0.06136, 0.05427, 0.05462, 0.0257, 0.00751, 0.00424, 0.00147, 0.0026];
    else
        N(m)=14;
        x(m,1:N(m)) = [0.34503, 0.44414, 0.57172, 0.72714, 0.88671, 1.27186, 1.78096, 2.43457, 3.47112, 4.4414, 6.14486, 8.70859, 12.04869, 17.28224];
        y(m,1:N(m)) = [0.70435, 0.23513, 0.23913, 0.23837, 0.24489, 0.23233, 0.21563, 0.16968, 0.104, 0.0711, 0.03218, 0.01289, 0.00323, 6.45752E-4];
        y_err(m,1:N(m)) = [0.28337, 0.03855, 0.05773, 0.05338, 0.06136, 0.05427, 0.05462, 0.0257, 0.00751, 0.00424, 0.00147, 0.0026, 0.00132, 2.69796E-4];
    end

    %Shao et al. (2011)
    m = 4; %sixth data set
    N(m) = 6;
    x(m,1:N(m)) = [0.71833, 1.09727, 1.67332, 2.64953, 4.55842, 7.04761];
    y(m,1:N(m)) = [0.3504, 0.17805, 0.27309, 0.24174, 0.06456, 0.01351];
    y_err(m,1:N(m)) = [0.02708, 0.01243, 0.01465, 0.00983, 9.11666E-4, 6.58086E-4];
    
    %Rosenberg et al. (2014)
    m = 5; %seventh data set
    if (PM10_only)
        N(m)=15;
        x(m,1:N(m)) = [0.55895, 0.67803, 0.85116, 1.05478, 1.34538, 1.7236, 1.98261, 2.14304, 2.35708, 2.57697, 2.80707, 3.02622, 3.32616, 4.67509, 7.66112];
        y(m,1:N(m)) = [1.60377, 1.37279, 0.95565, 0.85615, 0.22977, 0.26514, 0.1726, 0.36576, 0.29398, 0.25863, 0.18277, 0.20832, 0.08092, 0.02761, 0.01355];
        y_err(m,1:N(m)) = [0.42458, 0.27216, 0.18975, 0.09675, 0.03691, 0.04412, 0.02159, 0.06889, 0.01986, 0.04394, 0.02466, 0.08412, 0.02605, 0.00732, 0.00421];
    else
        N(m)=17;
        x(m,1:N(m)) = [0.55895, 0.67803, 0.85116, 1.05478, 1.34538, 1.7236, 1.98261, 2.14304, 2.35708, 2.57697, 2.80707, 3.02622, 3.32616, 4.67509, 7.66112, 11.6678, 16.7679];
        y(m,1:N(m)) = [1.60377, 1.37279, 0.95565, 0.85615, 0.22977, 0.26514, 0.1726, 0.36576, 0.29398, 0.25863, 0.18277, 0.20832, 0.08092, 0.02761, 0.01355, 0.00694, 0.00228];
        y_err(m,1:N(m)) = [0.42458, 0.27216, 0.18975, 0.09675, 0.03691, 0.04412, 0.02159, 0.06889, 0.01986, 0.04394, 0.02466, 0.08412, 0.02605, 0.00732, 0.00421, 0.00226, 0.0013];
    end
        
no_data_sets = m;