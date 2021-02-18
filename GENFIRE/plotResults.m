function plotResults(true_angles,refined_angles)
figure
plot(true_angles(:,1),'o')
hold on
plot(refined_angles(:,1),'r.')
legend('Original phi','Refined phi')

figure
plot(true_angles(:,2),'o')
hold on
plot(refined_angles(:,2),'r.')
legend('Original theta','Refined theta')

figure
plot(true_angles(:,3),'o')
hold on
plot(refined_angles(:,3),'r.')
legend('Original psi','Refined psi')
end