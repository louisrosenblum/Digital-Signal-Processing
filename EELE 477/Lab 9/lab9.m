%% Anthony Louis Rosenblum
%  Lab 9


% Anthony Louis Rosenblum
% EELE 477
% 3/7/2019
% Lab 9

%% Section 2-1

con2dis

%% Section 2-2

dconvdemo

%% Section 2-5

load echart.mat

% a

bdiffh = [1,-1]
yy1 = conv2(echart,bdiffh)

show_img(echart)
show_img(yy1)

%% b

yy2 = conv2(echart,bdiffh')

show_img(echart)
show_img(yy2)
%%

x1 = 256 * (rem(0:100,50)<10)

y1 = firfilt([1 -0.9],x1)

n = 0:100;

subplot(2,1,1)

stem(n,x1),xlim([0 75]),title("Input Signal"),ylabel("Magnitude"),xlabel("n")

subplot(2,1,2)

n2 = 0:101;

stem(n2,y1),xlim([0 75]),title("Output Signal"),xlabel("n"),ylabel("Magnitude")


%%
l = 0:22;
y2 = firfilt(0.9.^l,y1)
n3 = 0:123;
subplot(2,1,1)
stem(n2,y1),xlim([0 75]),title("W[n]"),ylabel("Magnitude"),xlabel("n")

subplot(2,1,2)
stem(n3,y2),xlim([0 75]),title("Y[n]"),ylabel("Magnitude"),xlabel("n")

%% 

xx1 = x1(1:51);
xx2 = y2(1:51);

yy1 = xx2 - xx1;

n = 0:50;

plot(n,yy1),title("Deconvolution Error"),ylabel("Magnitude"),xlabel("n")

%% 3-1-3

fs = 8000;

load labdat.mat

n1 = zeros(1601);
n(1) = 1;
n(1601) = 0.9;

y3 = firfilt(n1,x2)

soundsc(y3)

%% 3-2

x4 = zeros(100);
x4(1) = 1;

n = 0:99;
stem(n,x4)

w4 = firfilt([1 -0.9],x4);

l = 0:22;

r = 0.9
b = r.^l;

y4 = firfilt(b,w4)

n = 0:122

stem(n,y4),xlim([0 30]),ylim([-2 2]),title("Cascaded Output"),ylabel("Magnitude"),xlabel("n")

%% 3-2-2

load echart.mat

y1 = firfilt([1 -0.9],echart);

y2 = rot90(y1);

y3 = firfilt([1 -0.9],y2);

y4 = rot90(rot90(rot90(y3)));

show_img(y4)

l = 0:33;

r = 0.9
b = r.^l;

y5 = firfilt(b,y4)

y6 = rot90(y5);

y7 = firfilt(b,y6);

y8 = rot90(rot90(rot90(y7)));

show_img(y8)


