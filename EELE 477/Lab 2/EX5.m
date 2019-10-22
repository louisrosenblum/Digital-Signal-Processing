fk = [pi,pi,pi];
Xk = [2,2*exp(-1j * pi * 1.25),1-1j];
tstart = -0.5;
dur = (3 / pi);

[xx,tt] = syn_sin(fk, Xk, 10000, dur, tstart);
plot(tt,xx), xlabel('Time (s)'), ylabel('Value'),title('Figure 5-b'), text(-0.3302,1.639,'A = 1.639 Tm = 0.-0.3302'), text(0.3065,1.639, 'A =1.639 Tm = 0.3065')