%=====================================================================

% Purpose: To return the probability of detection
% when the target is located at a particular
% distance.
%
% Arguments
% Input(s):
% d - scalar or vector of distances from target to sensor(s)
% M - number of samples of recieved signal
% beta - threshold normalized over noise variance
% snr - signal to noise ratio when target and sensor are
% at a distance of 1 unit
% Output(s):
% result - Probability of detection for each element of the
% input distance vector 'd'
%
% Assumptions:
%
% Function Declaration:
% function result = P_d (d, M, beta, snr)
%=====================================================================
function result = P_d (d, M, beta, snr)
arg1 = 1 + snr ./ d.^2;
arg1 = beta ./ ( 2 * arg1 );
arg2 = M / 2;
result = 1 - gammainc (arg1,arg2);
%plot(d,result,'b*')
% ==============================================================
% Unused code
% ==============================================================
%
% one line implementation of Pd
% Pd(i) = 1 - gammainc( beta/( 2 * (1 + snr/d(i)^2) ), M/2 );k
% ==============================================================
%
% Implementation of Pd given in paper.
% Do not work as the gamma and incomplete gamma function given
% in paper and matlab are not compatible.
% Pd_den = gamma( M/2 );
% if Pd_den == 0
% error('Denomenator is zero')
% end
% arg1 = M/2;
% arg2 = 1 + snr/d(i)^2;
% arg2 = beta/2 * 1/arg2;
% Pd_num = gammainc( arg1,arg2 );
% Pd(i) = Pd_num / Pd_den;
% ==============================================================
%
% Non-matrix impelentation of the code (code with loops instead
% of matrix algebra.
% d_len = length(d);
% for i = 1:d_len
% arg1 = 1 + snr / d(i)^2;
% arg1 = beta / ( 2 * arg1 );
% arg2 = M / 2;
% result(i) = 1 - gammainc (arg1,arg2);
% end
