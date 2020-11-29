v = linspace(0,20,100);
vec_1=0.7*(1+zeros(1,100));
vec_2 = 0.6*(1+zeros(1,100));

alpha = vec_1./(1+exp(7*(v-15.25)));
mu = vec_2./(1+exp(7*(v-15.25)));
nu = vec_2./((1+exp((-7)*(v-12.75))).*(1+exp(7*(v-15.25))));

figure(1);
plot(v, alpha,'LineWidth',2); xlabel('v');ylabel('alpha');
legend ('alpha = 0.7/(1+exp(7*(v-15.25)))','Location','southwest');

figure(2);
plot(v,mu,'r','LineWidth',2);xlabel('v');ylabel('mu');
legend ('mu = 0.6/(1+exp(7*(v-15.25)))','Location','southwest');

figure(3);
plot(v,nu,'c','LineWidth',2);xlabel('v');ylabel('nu');
legend ('nu = 0.6//(1+exp(-7*(v-12.75)))*(1+exp(7*(v-15.25)))','Location','southwest');


