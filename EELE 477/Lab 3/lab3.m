%% Section 3 - Prelab

% a
z1 = 2 * exp(j*pi/3);
z2 = - sqrt(2) + 5*j;

zvect(z1), hold on, zcoords, ucplot, hold off;
zvect(z2), hold on, zcoords, ucplot, hold off;

zprint(z1);
zprint(z2);

%% b

z1con = conj(z1);
zprint(z1con);
z2con = conj(z2);
zprint(z2con)

z1in = 1/z1;
zprint(z1in)
z2in = 1/z2;
zprint(z2in)

%% c

zcat([1+j,-2+j,1-2*j])

%% d

z3 = z1 + z2;
zvect(z3), hold on, zcoords, ucplot;
zcat([z1,z2])
zprint(z3), hold off;

%% e

z4 = z1 * z2;
z5 = z2/z1;

zvect(z4), hold on, zcoords, ucplot;
zvect(z5), hold off;

zprint(z4)
zprint(z5)

%% f

subplot(2,2,1), zvect(z1), hold on, zcoords, ucplot;
zvect(z2);
zvect(z1+z2);

subplot(2,2,2), zvect(z2), hold on, zcoords, ucplot;
zvect(z2con);

subplot (2,2,3), zvect(z1), hold on, zcoords, ucplot;
zvect(z1in);

subplot (2,2,4), zvect(z1*z2), hold on, zcoords, ucplot;

%% 3-2

zdrill

%% 3-3

n = 200;
k = 1:n;
x = k/50;
rk = sqrt(x.*x+2.25)
sig = exp(j*2*pi*rk);
plot(x,real(sig))

%% Section 4

%4-1

[x1,t1] = one_cos(10^4,3*pi*10^6,-pi/4,10^-6)

plot(t1,x1),text(8.357e-8,1e4, "Peak")

%% 4-2

[xx0,tt0] = syn_sin([0,100,250],[10,14*exp(-pi/3),8*j],10000,0.1,0);

plot(tt0, real(xx0))

%% Section 5

[xx1,tt1] = syn_sin([50*pi,50*pi,50*pi],[-2,-1*exp(50*pi*-0.02*j),(2-j*3)],50*pi*2,3*(1/(50*pi)))

hold off;

plot(tt1,real(xx1)),text(0.007599,3.161,"Amplitude 3.161, peak at 0.007964"),
text(0.01398,3.161,"Peak at 0.00.01398"),text(0.001277,3.161,"Peak at 0.001277"),
title("Figure 5-a"),xlabel("Time (s)"), ylabel("Magnitude"), hold off

hand = (1-3j) * exp(2*pi*50*pi*j*tt1);
plot(tt1,real(hand)),text(0.007599,3.161,"Amplitude 3.161, peak at 0.007964"),
text(0.01398,3.161,"Peak at 0.00.01398"),text(0.001277,3.161,"Peak at 0.001277"),
title("Figure 5-c: Phasor Addition in MATLAB"),xlabel("Time (s)"), ylabel("Magnitude"), hold on

plot(tt1,cos(50*pi*2*pi*tt1-1.249)), text(0.01082,-1,"Hand calculation phase")

%% Section 6

y = 100;
d = 0.4;



f = 400;
amp = 1000;

t1 = delay_one(100)
t2 = delay_two(100)

tt1 = 0:1/1000000:3/400;
x1 = amp*cos((f*2*pi.*(tt1-t1)));
x2 = amp*cos((f*2*pi.*(tt1-t2)));

subplot(2,1,1),plot(tt1,x1),title("Receiver 1"),xlabel("Time (s)"), ylabel("Magnitude"),text(0.004266,1000,"0.004266"), hold on
subplot(2,1,2),plot(tt1,x2),title("Receiver 2"),xlabel("Time (s)"), ylabel("Magnitude"),text(0.003417,1000,"0.003417")


%% 6-1 theta vector

xv = -400:1:500;

f = 400;
amp = 1000;

dot = direction_finding(xv)

plot(xv,dot),title("Angle between Receivers"),xlabel("Vehicle Distance (m)"),ylabel('Radians')

hold off;


%% Function definitions

function[theta] = direction_finding(xv)

    dist1 = (xv - 0).^2 + (100-0)^2;
    dist2 = (xv - 0.4).^2 + (100-0)^2;
    d1 = sqrt(dist1);
    d2 = sqrt(dist2);
    t1 = d1/(333 + 1/3);
    t2 = d2/(333 + 1/3);
    
    t3 = t2 - t1;
    theta = t3*400 * 2 *pi;

end


function[t1] = delay_one(xv)

    dist = (xv - 0)^2 + (100-0)^2;
    d = sqrt(dist);
    t1 = d/(333 + 1/3);

end

function[t2] = delay_two(xv)

    dist2 = (xv - 0.4)^2 + (100-0)^2;
    d2 = sqrt(dist2);
    t2 = d2/(333 + 1/3);

end


function [xx,tt] = one_cos(a,w,p,dur)

nyq = w*16;
samp = 1/nyq;
tt = -dur/2:samp:dur/2;

xx = a * cos(w.*tt + p)

end

function [xx,tt] = syn_sin(fk, xk, fs, dur, tstart)

if nargin<5, tstart=0, end %--default value is zero

tt = tstart:dur/(fs):dur+tstart;

xx = xk*exp(j*2*pi*fk'*tt)

end
