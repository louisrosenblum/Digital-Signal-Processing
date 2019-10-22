%% Anthony Rosenblum
% Lab 4
% 1/31/2019

%% Pre-lab

fsamp = 11025;
dt = 1/fsamp;
dur = 1.8;
tt = 0 : dt : dur;
psi = 2*pi*(100 + 200*tt + 500*tt.*tt);
xx = real( 7.7*exp(j*psi) );
soundsc( xx, fsamp );

%a

length(tt); % 19846 points
% duration = 1.8s

%b

% tn = tt
% A = 7.7
% u = 500;
% f0 = 200
% pitchfork = 100

% c

length(psi)
psi(1) 
psi(19846)

% minimum = 628.32
% max = 1.3069e4

%d

% the signals frequency is increasing

%% Warm-up

% 3-1

beatcon

%% 3-2

[xx0,tt1] = mychirp(2500,500,1.5)

xx1 = real(xx0);

plot(tt1,xx1)

soundsc(xx1)


%% Section 4-1

[xx1,tt1] = beat(10,10,1000,10,11025,1)

plot(tt1,xx1),xlim([0 0.2]),title("Figure 4-1b"),xlabel("Time (s)"),ylabel("Amplitude");
%% Section 4-2

[xx2,tt2] = beat(10,10,2000,32,11025,0.26)

plot(tt2,real(xx2)),xlim([0 0.26]),title("Figure 4-2a"),xlabel("Time (s)"),ylabel("Amplitude");
%%


spectrogram(xx2,2048,[],2048,11025); colormap(1-gray(256)),title("Figure 4-2b"),xlabel("Frequency (hz)"),ylabel("Time (ms)");


%spectrogram(xx2,16,[],16,11025); colormap(1-gray(256)),title("Figure 4-2c"),xlabel("Frequency (hz)"),ylabel("Time (ms)");


%% Section 4-3

[xx3,tt3] = mychirp(4186.01,440,2.5,11025);

%plot(tt3,xx3)

soundsc(xx3)


spectrogram(xx3,1600,[],1600,11025); colormap(1-gray(256)),title("Figure 4-3"),xlabel("Frequency (hz)"),ylabel("Time (ms)");


%% Section 4-4

[xx4,tt4] = mychirp(3000,-2000,3,11025);

%plot(tt4,xx4),title("Figure 4-4"),xlabel("Time (s)"),ylabel("Amplitude");
soundsc(xx4)

spectrogram(xx4,2048,[],2048,11025); colormap(1-gray(256)),title("Figure 4-4 Spectrogram"),xlabel("Frequency (hz)"),ylabel("Time (ms)");

%% function defs

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

tt = tt1

end


function [xx,tt] = mychirp( f1, f2, dur, fsamp )
%MYCHIRP generate a linear-FM chirp signal
%
% usage: xx = mychirp( f1, f2, dur, fsamp )
%
% f1 = starting frequency
% f2 = ending frequency
% dur = total time duration
% fsamp = sampling frequency (OPTIONAL: default is 11025)
% xx = (vector of) samples of the chirp signal
% tt = vector of time instants for t=0 to t=dur
%

if( nargin < 4 ) %-- Allow optional input argument
fsamp =  11025;
end

tt = 0:dur/(fsamp*dur):dur;

w = f1 +(f2-f1)/dur .*tt;

plot(tt,w)

psi = f1*tt +(f2-f1)/(dur*2).*tt.*tt;

psirad = psi * 2 * pi;

xx = 20*cos(psirad);

end




