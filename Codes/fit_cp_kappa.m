cp = load('cp.txt');
kappa = load('kappa.txt');

figure(1)
plot(cp(:,2),cp(:,1),'bo','MarkerfaceColor','b');
xlabel('$\mathrm{Temperature~(K)}$','interpreter','latex');
ylabel('$\mathrm{Specific~Heat~(C_p,~J/kg/K)}$','interpreter','latex');

figure(2)
plot(kappa(:,2),kappa(:,1),'bo','MarkerfaceColor','b');
xlabel('$\mathrm{Temperature~(K)}$','interpreter','latex');
ylabel('$\mathrm{Thermal~Conductivity~(\kappa,~W/m/K)}$','interpreter','latex');


