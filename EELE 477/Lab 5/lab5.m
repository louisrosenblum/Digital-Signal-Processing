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

% 3-1

xx = key2note(10,40,3);

soundsc(xx)

%% 3-2


%--- play_scale.m
%---
scale.keys = [ 40 42 44 45 47 49 51 52 ];
%--- NOTES: C D E F G A B C
% key #40 is middle-C

scale.durations = 0.25 * ones(1,length(scale.keys));
fs = 11025; %-- or 8000 Hz
xx = zeros(1, sum(scale.durations)*fs+length(scale.keys) );
n1 = 1;
for kk = 1:length(scale.keys)
keynum = scale.keys(kk);
tone = key2note(10,scale.keys(kk),0.25)
n2 = n1 + length(tone) - 1;
xx(n1:n2) = xx(n1:n2) + tone; %<=== Insert the note
n1 = n2 + 1;
end
soundsc( xx, fs )

% 3-3

spectrogram(xx,512,[],512,fs,'yaxis')


%%

load bach_fugue.mat

key = theVoices.noteNumbers;

pulse = theVoices.startPulses;

dur = theVoices.durations;

pulse(2)

dur(2)

%% Section 4

% 4-3

load bach_fugue.mat

key = theVoices(1).noteNumbers;
pulse = theVoices(1).startPulses;
dur = theVoices(1).durations;

dur = dur/8;

fs = 11025;
xx = zeros(1,28*fs);
n1 = 1;
n2 = 0;

for kk = 1:length(key)
n1 = (pulse(kk)-1)/9*fs;
tone = key2note(10,key(kk),dur(kk))
tone = tone + key2note(10,key(kk)-12,dur(kk))
tone = tone + key2note(10,key(kk)+12,dur(kk))
tone = tone + key2note(10,key(kk)+24,dur(kk))
tone = tone + key2note(10,key(kk)-24,dur(kk))
tone = tone + key2note(10,key(kk)-36,dur(kk))
tone = tone + key2note(10,key(kk)+36,dur(kk))
n2 = n1 + length(tone) - 1;
xx(n1:n2) = xx(n1:n2) + tone; %<=== Insert the note
end

key = theVoices(2).noteNumbers;
pulse = theVoices(2).startPulses;
dur = theVoices(2).durations;
dur = dur/8;

n1 = 1;
n2 = 0;

for kk = 1:length(key)
n1 = (pulse(kk)-1)/9*fs;
tone = key2note(10,key(kk),dur(kk))
tone = tone + key2note(10,key(kk)-12,dur(kk))
tone = tone + key2note(10,key(kk)+12,dur(kk))
tone = tone + key2note(10,key(kk)+24,dur(kk))
tone = tone + key2note(10,key(kk)-24,dur(kk))
tone = tone + key2note(10,key(kk)-36,dur(kk))
tone = tone + key2note(10,key(kk)+36,dur(kk))
n2 = n1 + length(tone) - 1;
xx(n1:n2) = xx(n1:n2) + tone; %<=== Insert the note
end

key = theVoices(3).noteNumbers;
pulse = theVoices(3).startPulses;
dur = theVoices(3).durations;
dur = dur/8;

n1 = 1;
n2 = 0;

for kk = 1:length(key)
n1 = (pulse(kk)-1)/9*fs;
tone = key2note(10,key(kk),dur(kk))
tone = tone + key2note(10,key(kk)-12,dur(kk))
tone = tone + key2note(10,key(kk)+12,dur(kk))
tone = tone + key2note(10,key(kk)+24,dur(kk))
tone = tone + key2note(10,key(kk)-24,dur(kk))
tone = tone + key2note(10,key(kk)-36,dur(kk))
tone = tone + key2note(10,key(kk)+36,dur(kk))
n2 = n1 + length(tone) - 1;
xx(n1:n2) = xx(n1:n2) + tone; %<=== Insert the note
end

soundsc( xx, fs )

spectrogram(xx,3000,[],3000,fs,'yaxis'),title("Well-Tempered Clavier")

%%

soundsc( xx, fs )
%% function defs

function xx = key2note(X, keynum, dur)
% KEY2NOTE Produce a sinusoidal waveform corresponding to a
% given piano key number
%
% usage: xx = key2note (X, keynum, dur)
%
% xx = the output sinusoidal waveform
% X = complex amplitude for the sinusoid, X = A*exp(j*phi).
% keynum = the piano keyboard number of the desired note
% dur = the duration (in seconds) of the output note
%
fs = 11025; %-- or use 8000 Hz
tt = 0:(1/fs):dur;
num = 49 - keynum;
freq = 440 / 2^(1/12*num)
xx = real( X*exp(j*2*pi*freq*tt) );
end

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

