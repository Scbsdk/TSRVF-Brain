function q = f_to_srsf(f,time)
% f_to_srsf ---- Convert function to Square-Root Slope Function
%
% Usage: q = f_to_srsf(f,time)
%
% Input:
%   f: matrix of functions
%   time: column vector of time samples
%
% Output:
%   q: matrix of SRSFs


binsize = mean(diff(time));
fy = gradient(f,binsize);
q = fy./sqrt(abs(fy)+eps);

