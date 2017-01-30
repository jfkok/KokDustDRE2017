function [f] = dustPSD_power_law_fit (A, X_data)

%function that fits the dust PSD data to a power law following Kok (PNAS,
%2011)

f = A-2*X_data;
%diff_sq = (Y_data-calc).^2;
%chi_sq = sum(diff_sq)/size(X_data,2);
