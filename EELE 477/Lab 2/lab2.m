% Anthony Rosenblum
% EELE 477
% Lab 2
% 1/17/2019

%% Pre-lab, 3-1
j = sqrt(-1);

z1 = 10 * exp(j * 2/3 * pi);
z2 = -5 + 5*j;

zcat([j,01,-2j,1]);

z3 = z1 + z2;

zvect(z3), hold on, zcoords, ucplot, hold off;

z4 = z1 * z2;

zvect(z4), hold on, zcoords, ucplot, hold off;

z5 = z2/z1;

zvect(z5), hold on, zcoords, ucplot, hold off;

z6 = conj(z1);
z7 = conj(z2);

zvect(z6), hold on, zcoords, ucplot, hold off;
zprint(z6)

zvect(z7), hold on, zcoords, ucplot, hold off;
zprint(z7)

z8 = 1/z1;
z9 = 1/z2;

zvect(z8), hold on, zcoords, ucplot, hold off;
zprint(z8)

zvect(z9), hold on, zcoords, ucplot, hold off;
zprint(z9)

subplot(4,2,1), zvect(z1), hold, zvect(z2), hold on, zcoords, ucplot, hold off;

subplot(4,2,2), zvect(z6), hold, zvect(z7), hold on, zcoords, ucplot, hold off;

subplot(4,2,3), zvect(z8), hold, zvect(z9), hold on, zcoords, ucplot, hold off;

subplot(4,2,4), zvect(z4), hold on, zcoords, ucplot, hold off;

%% 3-2

zdrill

%% 3-3

%--- make a plot of a weird signal
%N = 200;
%for k=1:N
%xk(k) = k/50;
%rk(k) = sqrt( xk(k)*xk(k) + 2.25 );
%sig(k) = exp(j*2*pi*rk(k));
%end
%plot( xk, real(sig), "mo-" )


rk1 = sqrt(((1:N)/50).*((1:N)/50)+2.25) 
sig1 = exp(j*2*pi*rk1);
plot((1:N)/50, real(sig1), "mo")

%% 3-4

% Copied to function defs
%% 4-1

one_cos(95,200*pi,pi/5, 0.025)

%% 4-2

[xx0,tt0] = syn_sin([0,100,250],[10,14*exp(-j*pi/3),8*j],10000,0.1,0);

%% Section 5

fk = [pi,pi,pi];
Xk = [2,2*exp(-1j * pi * 1.25),1-1j];
tstart = -0.5;
dur = (3 / pi);

[xx,tt] = syn_sin(fk, Xk, 10000, dur, tstart);
plot(tt,xx), xlabel('Time (s)'), ylabel('Value'),title('Figure 5-b'), text(-0.3302,1.639,'A = 1.639 Tm = 0.-0.3302'), text(0.3065,1.639, 'A =1.639 Tm = 0.3065'), text(-0.0131,1.639,'A=1.639 Tm=-0.0131')

%% 5-c
t = -0.5:1000:3/pi-0.5;
x1 = 2*exp(j*pi*t);
x2 = 2*exp(-1j * pi * 1.25) * exp(j*pi*t);
x3 = (1-1j)* exp(j*pi*t);

x4 = real(x1) + real(x2) +  real(x3)
x5 = imag(x1) + imag(x2) +  imag(x3)

amp = sqrt((x4)^2 + (x5)^2)
phase = atan(-x4/x5)

%% Section 6

yt = 1500;
bx = 100;
by = 900;

t1 = direct(yt,0)

t2 = bounce(yt,bx,by,0)

f = 200 * 1000000;
per = 1/f;
w = f * 2 * pi;
tt = -2*per:1/100000000000:2*per;
r = cos(w*(tt-t1)) + cos(w*(tt-t2));

plot(tt,r), xlabel('Time (s)'), ylabel('Value'),title('Figure 6-c'),text(5.3e-10,1.581,'A = 1.581 Tm = 5.3e-10')
%r = exp(-t1/per*2*pi*j)*exp(w*tt*j) + exp(-t2/per*2*pi*j)*exp(w*tt*j);

%plot(tt,real(r)),xlabel('Time (s)'), ylabel('Value'),title('Figure 6-d'), text(5.3e-10,1.581,'A = 1.581 Tm = 5.3e-10')

%% 6-e

xv = 0:300;
xx3 = go(xv);

plot(xv,real(xx3)),xlabel('Vehicle X Distance (m)'), ylabel('Strength'), title('Figure 6-f')

%% Function defs

function[xx3,t1,t2] = go(xv)
    t1 = sqrt((2000)^2 + (xv).^2) / (3 * 10^8);
    step = sqrt((2000-1000)^2 + (200)^2) + sqrt((200-xv).^2 + (1000)^2);
    t2 = step / (3 * 10^8);

    f = 150 * 1000000;
    per = 1/f;
    w = f * 2 * pi;
    tt = -10*per:1/100000000000:10*per;
    
    
    xx3 = exp(-t1/per*2*pi*j).*exp(w*xv*j) + exp(-t2/per*2*pi*j).*exp(w*xv*j);
end

function[xx]= direct(yt,xv)
    xx = sqrt((yt)^2 + (xv)^2) / (3 * 10^8);

end

function[xx2] = bounce(dt,bx,by,xv)
    step = sqrt((dt-by)^2 + (bx)^2) + sqrt((bx-xv)^2 + (by)^2);
    xx2 = step / (3 * 10^8);
end





function [xx,tt] = one_cos(amp, frq, phs, dur)
tt = 0:1/(100*frq):dur;
xx = amp * cos(frq*tt+phs);

plot(tt,xx)
end

function [xx,tt] = goodcos(ff,dur)
tt = 0:1/(100*ff):dur; %-- gives 100 samples per period
xx = cos(2*pi*ff*tt);

end

function [xx,tt] = syn_sin(fk, Xk, fs, dur, tstart)
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

xx = real(Xk * exp(p))

%subplot(2,1,1),plot(tt,real(xx))
%subplot(2,1,2),plot(tt,angle(xx))


end







