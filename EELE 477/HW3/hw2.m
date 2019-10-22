tt = 0:1/100:3;

x = 3 * sin(1/6*2*pi*tt);

y = 0.5 * sin(8*2*pi*tt);

z = x + y;

plot(tt,z)

%%

t1 = 2:1/100:5;

x1 = 4 * exp(log(1/4)*(t1-3));

plot(t1,x1)