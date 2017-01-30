function [intercept, slope, SE_intercept, SE_slope, covariance] = linearfit (x, y, y_err)

%this function performs a linear fit to the x, y, and y_err data

N = size(x,2); %total number of measurements

if (y_err==0)
    known_uncertainty = false; %uncertainties of the measurements are not known. they will therefore be assumed to all be equal, with a value extracted from the least squares fit.
else
    known_uncertainty = true; %uncertainties of the measurements are known and will be used in the fitting procedure
end

if (known_uncertainty)
    delta = sum(1./y_err.^2)*sum(x.^2./y_err.^2)-sum(x./y_err.^2).^2;
    intercept = (1/delta)*(sum(x.^2./y_err.^2)*sum(y./y_err.^2)-sum(x./y_err.^2)*sum(x.*y./y_err.^2)); %intercept = a in Bevington and Robinson
    slope = (1/delta)*(sum(1./y_err.^2)*sum(x.*y./y_err.^2)-sum(x./y_err.^2)*sum(y./y_err.^2));
    SE_intercept = sqrt((1/delta)*sum(x.^2./y_err.^2));
    SE_slope = sqrt((1/delta)*sum(1./y_err.^2));
    %covariance = sqrt((1/delta)*sum(x./y_err.^2));
    for p=1:N
        del_int(p) = (1/delta)*(1/y_err(p)^2)*(sum(x.^2./y_err.^2)-x(p)*sum(x./y_err.^2)); %intercept = a in Bevington and Robinson
        del_slope(p) = (1/delta)*(1/y_err(p)^2)*(x(p)*sum(1./y_err.^2)-sum(x./y_err.^2)); %slope = b in Bevington and Robinson
    end
    SE_intercept_alt = sqrt(sum(y_err.^2.*del_int.^2));
    SE_slope_alt = sqrt(sum(y_err.^2.*del_slope.^2));
    covariance = sqrt(sum(y_err.^2.*del_slope.*del_int));
else
    delta = N*sum(x.^2)-(sum(x))^2;
    intercept = (1/delta)*(sum(x.^2)*sum(y)-sum(x)*sum(x.*y)); %intercept = a in Bevington and Robinson
    slope = (1/delta)*(N*sum(x.*y)-sum(x)*sum(y)); %slope = b in Bevington and Robinson
    y_err = sqrt((1/(N-2))*sum((y-intercept-slope*x).^2));
    SE_intercept = sqrt(y_err^2*sum(x.^2)/delta);
    SE_slope = sqrt(N*y_err^2/delta);
    for p=1:N
        del_int(p) = (1/delta)*(sum(x.^2)-x(p)*sum(x)); %intercept = a in Bevington and Robinson
        del_slope(p) = (1/delta)*(N*x(p)-sum(x)); %slope = b in Bevington and Robinson
    end
    SE_intercept_alt = sqrt(sum(y_err.^2.*del_int.^2));
    SE_slope_alt = sqrt(sum(y_err.^2.*del_slope.^2));
    covariance = sqrt(sum(y_err.^2.*del_slope.*del_int));
end

1;

%for i=1:size(intercept,2)
%    for j=1:size(slope,2)
        %for p=1:size(x,2)
            %del_int(p) = (1/delta)*(1/y_err(p)^2)*(sum(x.^2./y_err.^2)-x(p)*sum(x./y_err(p).^2)); %intercept = a in Bevington and Robinson
            %del_slope(p) = (1/delta)*(1/y_err(p)^2)*(x(p)*sum(1./y_err.^2)-sum(x./y_err.^2)); %slope = b in Bevington and Robinson
            %calc(i,j,p) = intercept(i)+slope(j)*x(p);
            %diff(i,j,p)=(1/y_err(p))^2*(calc(i,j,p)-y(p))^2; %change x#
        %end
%        chi_square(i,j)=sqrt(sum(diff(i,j,:)));
%    end
%end




%chi_square;
%figure(2); clf; contour(intercept,slope,chi_square);
%[a,j_final]=min(min(chi_square));
%[a,i_final]=min(chi_square(:,j_final));
%chi_square(i_final, j_final);

%intercept_alt = intercept(i_final)
%slope_alt = slope(j_final)
