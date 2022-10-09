%������չϵ��C_fcg������������㸯ʴƣ������
clear; close all; clc; format compact; format long;
cd 'E:\CYW\Matlab\CF';      % Path of Codes

C_fcg = randn(10,1)*2.605e-13+5.21e-13;        % FCG rate in air   ƽ��ֵ��5.21e-13������ϵ����50%
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


