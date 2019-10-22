function [xx,tt] = syn_sin(fk, xk, fs, dur, tstart)

if nargin<5, tstart=0, end %--default value is zero

tt = tstart:dur/(fs):dur+tstart;

xx = xk*exp(j*2*pi*fk'*tt)

end


function [xx, tt] = beat(A, B, fc, delf, fsamp, dur)
% A = amplitude of lower frequency cosine
% B = amplitude of higher frequency cosine
% fc = center frequency
% delf = frequency difference
% fsamp = sampling rate
% dur = total time duration in seconds
% xx = output vector of samples
%--Second Output:
% tt = time vector corresponding to xx

tt = 0:1/fsamp:dur;

fk = [fc - delf, fc + delf];
xk = [A,B]

[xx,tt1] = syn_sin(fk,xk,fsamp,dur,0);

end