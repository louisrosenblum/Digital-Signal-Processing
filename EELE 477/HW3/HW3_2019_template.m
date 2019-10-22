%clear all
%close all

To = 5;
%--------------------------------
% Segment 1 
%--------------------------------
t1=0:1/100:3;
y11=3*sin(2*pi*(1/6)*t1);
y12=0.5*sin(2*pi*8*t1);
y13 = y11+y12;
h=plot(t1,y13,'g'); hold on
set(h,'LineWidth',3.0)
h=plot(t1-To,y13,'g'); hold on
set(h,'LineWidth',3.0)
h=plot(t1+To,y13,'g'); hold on
set(h,'LineWidth',3.0)
plot(t1,y11,'g:')
plot(t1-To,y11,'g:')
plot(t1+To,y11,'g:')

%--------------------------------
% Segment 2
%--------------------------------
t2 = 3:1/100:4;
a = log(1)-log(4);
b = log(4) - 3*a;
y2 = exp(a*t2+b);
h=plot(t2,y2,'g');
set(h,'LineWidth',3.0)
h=plot(t2-To,y2,'g');
set(h,'LineWidth',3.0)
h=plot(t2+To,y2,'g');
set(h,'LineWidth',3.0)

%--------------------------------
% Segment 3
%--------------------------------
t4 = 4:1/100:5;
y4 = ones(size(t4))*3;
h=plot(t4,y4,'g');
set(h,'LineWidth',3.0)
h=plot(t4-To,y4,'g');
set(h,'LineWidth',3.0)
h=plot(t4+To,y4,'g');
set(h,'LineWidth',3.0)

%--------------------------------
% Dotted Lines
%--------------------------------
for kv=-1:1
    offset = kv*To;
    plot([0 5]+offset,[-1 -1],'g:')
    plot([0 5]+offset,[0 0],'g:')
    plot([0 5]+offset,[1 1],'g:')
    plot([0 5]+offset,[2 2],'g:')
    plot([0 5]+offset,[3 3],'g:')
    plot([0 5]+offset,[4 4],'g:')
    plot([1 1]+offset,[-2 5],'g:')
    plot([2 2]+offset,[-2 5],'g:')
    plot([3 3]+offset,[-2 5],'g:')
    plot([4 4]+offset,[-2 5],'g:')
    plot([3 3]+offset,[0 4],'g')
    plot([4 4]+offset,[1 3],'g')
    plot([5 5]+offset,[0 3],'g')
end

axis([-To 2*To -0.5 4.5])

%h=title('Waveform Template (T_o = 5 seconds)')
%set(h,'FontSize',18)
h=xlabel('Time (seconds)');
set(h,'FontSize',18)

set(gca,'FontSize',18)

% hw 3

%%

k = 1:1:1000;
fk = k/5;

z = pi * k/5;

a1 = (3 * cos(pi-6*z) - 3)./(2*z-pi/3); %a1 and a2 coefficients from hand calculations segment 1
a2 = (0.5 * cos(48*pi - 6*z) - 0.5)./(2*z - 16*pi);
afix = (0.5 * 6 * pi/5*sin(48*pi-6*40*pi/5))/(2*pi/5);
if length(a2) >= 40
    a2(40) = afix;
end

a3 = a1 + a2;

c1 = 3./(-2*j*z).*(exp(-5*j*2*z)-exp(-4*j*2*z)); %c1 hand calculations from segment 3

b1 = 4 * exp(log(1/4)*(-3))./(log(1/4) - 2*j*z).*(exp(4*(log(1/4)-2*j*z))-exp(3*(log(1/4) -2*j*z))); %b1 hand coefficients from hand calculations segment 2
b2 = b1/(2.342) + 1.2;
xk = 2/5*(a3+b1+c1);

[xx2,tt1] = sin_sin(fk,b1,1000,15,-5);

ack1 = -(225*exp(-(6*pi*j*k)/5).*(exp((6*pi*j*k)/5)+1))./(pi*(36*j^2*k.^2+25)); %ack1 and 2 verified coefficients on integralcalculator.com

ack2 = (50*exp(-(6*pi*j*k)/5).*(exp((6*pi*j*k)/5)-1))./(pi*(j^2*k.^2+1600));
if length(ack2) >= 40
    ack2(40) = afix;
end

tt2 = tt1;

tt2(1:3000) = 1;
tt2(3001:5000) = 0;
tt2(5001:8000) = 1;
tt2(8001:10000) = 0;
tt2(10001:13000) = 1;
tt2(13001:15000) = 0;

tt3 = tt1;

tt3(1:3000) = 0;
tt3(3001:4000) = 1;
tt3(4001:8000) = 0;
tt3(8001:9000) = 1;
tt3(9001:13000) = 0;
tt3(13001:14000) = 1;
tt3(14001:15000) = 0;

xx22 = xx2/2.79 + 0.58;

xx22 = xx22 .*tt3; % good
a1 = ack1;
a2 = ack2;


[xx3,tt1] = sin_sin(fk,c1,1000,15,-5);
xx32 = xx3/2.5 + 0.6; % good

[xx1,tt1] =  sin_sin(fk,a1,1000,15,-5);

xx12 = xx1/2.5 * -1 + 1.15;


[xx0,tt1] =  sin_sin(fk,a2,1000,15,-5);


xk = -1*a1 + a2 + b1 + c1;
xx01 = xx0 .* tt2;

xx4 = xx22 + xx32+xx12+xx01;


xx_plot = plot(fk,angle(xk)),set(xx_plot,'LineWidth',1.0),title("Frequency Spectrum: 1000 terms"),ylabel("Angle"),xlabel('Frequency (hz)')


%% Function definitions


function [xx,tt] = syn_sin(fk, xk, fs, dur, tstart)

if nargin<5, tstart=0, end %--default value is zero

tt = tstart:dur./(fs):dur+tstart;

xx = xk*exp(j*2*pi*fk'*tt);

end

function [xx,tt] = sin_sin(fk, Xk, fs, dur, tstart)
%SYN_SIN Function to synthesize a sum of cosine waves
% usage:
% [xx,tt] = syn_sin(fk, Xk, fs, dur, tstart)
% fk = vector of frequencies
% (these could be negative or positive)
% Xk = vector of complex amplitudes: Amp*eˆ(j*phase)
% fs = the number of samples per second for the time axis
% dur = total time duration of the signal
% tstart = starting time (default is zero, if you make this input optional)
% xx = vector of sinusoidal values
% tt = vector of times, for the time axis
%
% Note: fk and Xk must be the same length.
% Xk(1) corresponds to frequency fk(1),
% Xk(2) corresponds to frequency fk(2), etc.

if nargin<5, tstart=0, end %--default value is zero

tt = tstart:1/fs:(tstart+dur);

j = sqrt(-1);

fk = fk*2*pi;

p = fk' * tt *j;

xx = real(Xk * exp(p));

end

