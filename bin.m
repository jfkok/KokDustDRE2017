function [x, Px, x_neg2sd, x_neg1sd, x_median, x_pos1sd, x_pos2sd] = bin (var, bin_size, filetext, log_scale, process_log)

%function that sorts the data in var into bins of size bin_size, and writes
%out the result to filetext

write_data = false;

if (log_scale)
    var = log10(var(find(var>0)))+process_log; %process_log equals the number of orders of magnitude to be added to the variable such that all logs of the variable are > 0. For the reptation length, this is usually 9, such that the units on the reptation length are nm
    %min_bin = 1; %this worked for the creep length for some reason, but not for the impact speed...
    min_bin = round(min(min(var))/bin_size-0.5);
    max_bin = round(max(max(var))/bin_size+0.5);
else
    min_bin = round(min(min(var))/bin_size-0.5);
    max_bin = round(max(max(var))/bin_size+0.5);
end
bin_content = zeros(-min_bin+max_bin,1);
total_no = 0;

for r=1:size(var,1)
    for c=1:size(var,2)
        if (var(r,c)~=0)
            bin_no = -min_bin + round(var(r,c)/bin_size + 0.5);
            bin_content(bin_no) = bin_content(bin_no) + 1;
            total_no = total_no + 1;
        end
    end
end
x = ((1:size(bin_content,1))-0.5+min_bin)*bin_size;

if (log_scale) %for when the pdf is in log space
    bin_low = 10.^(x-bin_size/2); %the lower limit on the bin in linear coordinates
    bin_hi = 10.^(x+bin_size/2); %the upper limit on the bin in linear coordinates
    Px = bin_content'./(total_no*(bin_hi-bin_low)); %the probability density = number of occurrences / (total number of occurrences * width of the bin in linear space)
else %for when the pdf is in linear space
    Px = (bin_content/(total_no*bin_size))'; %the probability density = number of occurrences / (total number of occurrences * width of the bin in linear space)
end
sum_Px = sum(Px);
%figure(1); clf; plot(x-0.5*bin_size,Px);

for r=1:size(x,2)
    cum_Px(r) = sum(Px(1:r))/sum_Px;
    if (r==1)
        if (cum_Px(r) > 0.02275)
            x_neg2sd = x(r);
        end
        if (cum_Px(r) > 0.1587)
            x_neg1sd = x(r);
        end
        if (cum_Px(r) > 0.5)
            x_median = x(r);
        end
        if (cum_Px(r) > 0.8413)
            x_pos1sd = x(r);
        end
        if (cum_Px(r) > 0.97725)
            x_pos2sd = x(r);
            break;
        end
    else
        if (cum_Px(r) >= 0.02275 && cum_Px(r-1) < 0.02275)
            x_neg2sd = x(r);
        end
        if (cum_Px(r) >= 0.1587 && cum_Px(r-1) < 0.1587)
            x_neg1sd = x(r);
        end
        if (cum_Px(r) >= 0.5 && cum_Px(r-1) < 0.5)
            x_median = x(r);
        end
        if (cum_Px(r) >= 0.8413 && cum_Px(r-1) < 0.8413)
            x_pos1sd = x(r);
        end
        if (cum_Px(r) >= 0.97725 && cum_Px(r-1) < 0.97725)
            x_pos2sd = x(r);
        end
    end
end

if (log_scale)
    bin_width=10.^(x+bin_size/2)-10.^(x-bin_size/2); %the bin-width in linear coordinates
    Px=Px/sum(Px.*bin_width); %normalizing such that the integral over Px yields 1 in linear coordinates (note that Px is not unitless, so for converting from nm to m in creep length would require multiplication by 10^9 for Px
    x = 10.^x; %converting x back to linear space
end
%Px = Px/sum(Px);

if (write_data == true)
    fidData = fopen(strcat('bin',filetext,'.txt'), 'wt');
    %fprintf(fidData,'x      Px\n');
    for r = 1:size(bin_content,1)
        if (log_scale)
            fprintf(fidData, '%1.5e %1.5e %1.5e \n', x(r), bin_width(r), Px(r));
        else
            fprintf(fidData, '%1.5e %1.5e %1.5e \n', x(r), bin_size, Px(r));
        end
    end
    fclose(fidData);
end
            
            
        