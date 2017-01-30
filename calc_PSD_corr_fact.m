function [PSD_corr_fact, y_V_corr, y_err_V_corr, PSD_emit_avg] = calc_PSD_corr_fact(x, y_N, y_err_N, weight, D, PSD_lnV_emit_norm, j_end)

%this function calculates the factor by which each data sets needs to be
%corrected to maximize agreement against the fit, which is necessary
%because the non-continuity of dust PSD measurements in D make it such that
%there is no single normalization factor

%converting dN/dlnD measurements to dV/dlnD measurements
y_V = y_N.*(3.14159/6).*x.^3;
y_err_V = y_err_N.*(3.14159/6).*x.^3;

%obtaining the mean bootstrap PSD
for p = 1:size(D,2)
    A=sort(PSD_lnV_emit_norm(:,p));
    PSD_emit_avg(p) = A(round(0.5*j_end)); %mean bootstrap PSD
end

C = logspace(-3,1,100);

for i=1:size(C,2)
    for p=1:size(x,2)
        j = max(find(D<x(p)));
        diff(i,p) = weight(p)*(log(C(i)*y_V(p))-log(PSD_emit_avg(j)))^2;        
    end
    chi_square(i) = sum(diff(i,:));
end

[A,i_min] = min(chi_square);
C_min = C(i_min);
delta_C = C_min/200;
C = C(i_min-1):delta_C:C(i_min+1);

clear chi_square
for i=1:size(C,2)
    for p=1:size(x,2)
        j = max(find(D<x(p)));
        diff(i,p) = weight(p)*(C(i)*log(y_V(p))-log(PSD_emit_avg(j)))^2;
    end
    chi_square(i) = sum(diff(i,:));
end
[A,i_min] = min(chi_square);
PSD_corr_fact = C(i_min);

%obtaining the corrected data
y_V_corr = y_V*PSD_corr_fact;
y_err_V_corr = y_err_V*PSD_corr_fact;
1;
