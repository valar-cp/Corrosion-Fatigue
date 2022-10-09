%生成扩展系数C_fcg随机数，并计算腐蚀疲劳寿命
clear; close all; clc; format compact; format long;
cd 'E:\CYW\Matlab\CF';      % Path of Codes

C_fcg = randn(10,1)*2.605e-13+5.21e-13;        % FCG rate in air   平均值：5.21e-13，变异系数：50%
C_fcg = C_fcg(C_fcg>8.3e-14);
x = length(C_fcg);

for i=1:x
    t_sum = DemoTest(C_fcg(i));
    m(i,1) = t_sum;
end
figure(1);
histogram(C_fcg);
figure(2);
histogram(m);


