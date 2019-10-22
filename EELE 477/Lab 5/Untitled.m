% Anthony Rosenblum
% DSP - Lab 5
% 2/7/2019

%% Pre-lab

% 2-2a

tt = 0:0.5/11025:0.5;
x1 = 100 * cos(2*pi*800.*tt -pi/3)

soundsc(x1)
%% 2-2b
tt2 = 0:0.8/11025:0.8;
x2 = 80 * cos(2*pi*1200.*tt2+pi/4)

soundsc(x2)

%% 2-2c

xx = [x1,zeros(1,1103),x2]

soundsc(xx)



%% 2-2d
ttx = (1/11025)*(1:length(xx));
plot(ttx,xx)

%% 2-2e

soundsc(xx,22050)

% The duration of the signal is shortened, and the pitch is increase

%% 2-3

x.Amp = 7;
x.phase = -pi/2;
x.freq = 100;
x.fs = 11025
x.timeInterval = 0:(1/x.fs):0.05;
x.values = x.Amp*cos(2*pi*(x.freq)*(x.timeInterval) + x.phase);
x.name = "My Signal";
x %---- echo the contents of the structure "x"
plot( x.timeInterval, x.values )
title( x.name )

%% 2-4

dbstop if error
[xn,tn] = coscos(2,3,20,1)

%% 3 - Warm Up



%% function defs
function [xx,tt] = coscos( f1, f2, fs, dur )
% COSCOS multiply two sinusoids
%
t1 = 0:(1/fs):dur;
t2 = 0:(1/f2):dur;
cos1 = cos(2*pi*f1*t1);
cos2 = cos(2*pi*f2*t2);
xx = cos1.* cos2;
tt = t1;

end

function [xx,tt] = syn_sin(fk, xk, fs, dur, tstart)

if nargin<5, tstart=0, end %--default value is zero

tt = tstart:dur/(fs):dur+tstart;

xx = xk*exp(j*2*pi*fk'*tt)

end

