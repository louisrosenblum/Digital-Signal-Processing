% hw 3

k = 1:5;

xk = 2/5 .* ((3.*cos(pi-6.*pi.*k./5)-3)./(2.*pi.*k./5-pi./3) + (0.5.*cos(48.*pi-6.*pi.*k./5)-0.05)./(2.*pi.*k./5 - 16 .* pi) + 4.*(exp(log(1./4).*-3))./(log(1./4)-2.*pi.*j.*k./5).*(exp(4.*(log(1./4)-2.*pi.*j.*k./5)) -exp(3.*(log(1./4)-2.*pi.*j.*k./5))) + 3./(-j.*2.*pi.*k./5).*(exp(-5.*j.*2.*pi.*k./5) -exp(-4.*j.*2.*pi.*k./5)));




%% Function definitions


function [xx,tt] = syn_sin(fk, xk, fs, dur, tstart)

if nargin<5, tstart=0, end %--default value is zero

tt = tstart:dur./(fs):dur+tstart;

xx = xk*exp(j*2*pi*fk'*tt)

end
