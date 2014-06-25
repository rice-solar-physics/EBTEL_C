%adaptive_test_plot.m

clear all
close all

test_array = load('adaptive_testing.txt');

x = test_array(:,2);
t = test_array(:,1);
tau = test_array(:,3);

%Plot the data
figure
hold on
plot(t,x/max(x)*max(tau),'.')
plot(t,tau,'r--')
legend('x','tau','Location','Best')