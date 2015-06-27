function SEM=compute_standard_error(vals)

%COMPUTE_STANDARD_ERROR computes standard error of the mean

standard_deviation=nanstd(vals);
len_vals=length(vals);

SEM=standard_deviation./sqrt(len_vals);


end