
run Typ3ModellParameter;

close all;

sim('Typ3SimulinkModell'); 

figure;
subplot(3,1,1);
plot(tOutSwitching,dOutSwitching,'b',t_sim,d_Simulink,'r','LineWidth',2);
legend('d_{mCode}','d_{Simulink}');
xlabel('time [s]');
ylabel('duty cycle [-]');
title('Switching dc-to-dc converter');
grid;
subplot(3,1,2);
plot(tOutSwitching,xOutSwitching(1,:),'b',t_Simulink,iL1_Simulink,'r','LineWidth',2);
legend('i_{L1 mCode}','i_{L1 Simulink}');
xlabel('time [s]');
ylabel('current [A]');
grid;
subplot(3,1,3);
plot(tOutSwitching,xOutSwitching(4,:),'b',t_Simulink,vC_Simulink,'r','LineWidth',2);
legend('v_{C mCode}','v_{C Simulink}');
xlabel('time [s]');
ylabel('voltage [V]');
grid;