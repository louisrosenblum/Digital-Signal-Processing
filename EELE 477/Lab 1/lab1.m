% Anthony Rosenblum
% Lab 1
% EELE 477 - DSP

%% Pre-lab 

% 1.3e
pi * pi - 10
sin(pi/4)
ans ^ 2

% 1.3f
x = sin(pi/5);
cos(pi/5)
y = sqrt(1 - x*x)
ans

% 1.3g
z = 3 + 4i
w = -3 + 4j
real(z), imag(z)
abs([z,w])
conj(z+w)
angle(z)
exp(j*pi)
exp(j*[pi/4,0,-pi/4])

%% 2 - Warm Up

% 2.1a

jkl = 0:6
% Vector of ints from 0 to 6

jkl = 2:4:17
% Vector of ints spaced by 4 from 2 to 17 (2,6,10,14)

jkl = 99:-1:88
% Vector of ints counting down from 99 to 88

ttt = 2: (1/9): 4
% Vector from 2 to 4 of values spaced by 1/9

tpi = pi * [0:0.1:2]
% Vector from 0 to 2 of values spaced by 0.1, then multiplied by a constant
% of pi

%% 2.1b

xx = [zeros(1,3),linspace(0,1,5),ones(1,4)]
% Long vector of zeros, followed by linspace (5 values from 0 to 1 equally
% spaced), followed by ones

xx(4:6)
% indices of xx

size(xx)
% 1 by 12 matrix 

length(xx)
% vector of length 12

xx(2:2:length(xx))
% xx indices 2 to 12, every other (even indices only)

%% 2-1c
 yy = xx; 
 
 yy(4:6) = pi*(1:3);
 
 xx(2*(1:6)) = pi^pi

 %% 2-2a
 
 xk = cos(pi*(0:11)/4)
 % Cos of every value in vector is calculated and stored in xk
 
 xk(1)
 % xk(1) is 1
 
 % xk(0) is undefined, MATLAB starts indexing at 1????
 
 %% 2-2b
 
 yy = [];
 for k=-5:5
        yy(k+6) = cos(k*pi/3)
 end
 yy
 
 % It is necessary to write +6, because the indices must start at 1, not
 % -5. But we want to use the k value -5 in our evaluation.
 % Without the +6 you will get an indice error
 
 %% 2-1c
 
 x = [-3 -1 0 1 3];
 
 y = x.*x - 3*x;
 
 plot(x,y);
 
 z = x +y*sqrt(-1)
 plot(z)
 
 %% 2-d
 
 tt = -1:0.01:1;
 % time vector
 
 xx = cos(5*pi*tt);
 % cos function
 
 zz = 1.4*exp(j*pi/2)*exp(j*5*pi*tt);
 
 % amp mult by 1.4, phase shift of pi/2, multiplied by function
 % cos(5pi*tt) + jsin(5pi*tt)
 
 plot(tt,xx,'b-', tt, real(zz), 'r--'), grid on
 % Real zz is a sinusoid because the real part of the function zz is a cos
 % wave which was already shifted by pi/2 phase and amp multiplicaiton of
 % 1.4
 
 title('TEST PLOT of a SINUSOID')
 xlabel('TIME (sec)')

%% 2-3

t = 0:(1/11025):0.9;

m = sin(1600 * 2*pi*t)

soundsc(m)

%% 3 Lab Excercise

tt = -1/4000:1/(4000*100):1/4000;

k = (37.2/6)*1/4000;
p = -(41.3/10)*1/4000;

n = 19*cos(2*pi*4000*(tt-k));

m = 19*1.2*cos(2*pi*4000*(tt-p));

o = n + m;

p1 = k*4000*2*pi 
p2 = p * 4000*2*pi

subplot(3,1,1), plot(tt,n), title('Anthony Rosenblum - X1'), xlabel('Time (s)'), ylabel('Value'), text(-0.0002,19,'(-0.0002,19) A = 19 Tm = 0.00155 Phase: -5.026 radians')
subplot(3,1,2), plot(tt,m), title('X2'), xlabel('Time (s)'), ylabel('Value'), text(-0.0000325,22.76,'(-3.25e-5,22.76) A = 22.8 Tm = -0.0010325 Phase: -0.8169 radians')
subplot(3,1,3), plot(tt,o), title('X3'), xlabel('Time (s)'), ylabel('Value'), text(0,21.53,'(0,21.53) A = 21.53 Tm = 2.68e-6  Phase: 0.06747 radians')

%% 3-2

x = real(19*exp(sqrt(-1)*4000*2*pi*(tt))*exp(sqrt(-1)*(5.026)))

plot(tt,x), title('One Line of Code Complex Exponential'),xlabel('Time (s)'), ylabel('Value')












