close all;
clear all;

rng(123);

% Process Parameters

vel = 500; % scan speed (mm/s)
pow = 160; % power (W)
pre = 650; % Pre-heat temperature (C)

% Material Properties

S = 825; % yield strength (MPa)
E = 110; % elastic strength (GPa)
rho = 4428; % density (kg/m3)

C0_Cp = 540; C1_Cp = 0.43; C2_Cp = -0.000032;
C0_k = 7.2; C1_k = 0.011; C2_k = 0.0000014;

m = 12;
k = m+1;
alpha = 3;
M = floor( alpha *k * log(m) );

nom = [vel,pow,pre,S,E,rho,C0_Cp,C1_Cp,C2_Cp,C0_k,C1_k,C2_k];

L = zeros(1,m); U = zeros(1,m);
L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);

N = 20; % 20 samples in the 12D parameter space

% generate 20 samples in [-1,1]
xi  = -1 + 2*lhsdesign(N,m);
y = -1 + 2 * lhsdesign(M,m);

%pts_x = gen_samples(xi,L,U);
[lambda,W] = compute_subspace(m,N,M,xi,y);

function pts_x = gen_samples(xi,L,U)
pts_x = zeros(size(xi,1),size(xi,2));

for i = 1:size(pts_x,1)
  for j = 1:size(pts_x,2)
    pts_x(i,j) = L(1,j) + 0.5*(U(1,j)-L(1,j)).*(xi(i,j)+1);
  end
end

% Save physical points to a file
save('pts_gradfree_N20.txt','pts_x','-ASCII');

end 

function [lambda,W] = compute_subspace(m,N,M,xi,y)

Zf = load('Zf.txt');
Zf3 = Zf(:,3);

f = zeros(N,1);
f(:,1) = Zf3(1:N);

xi_ref = xi;
y_ref = y;

%
% find the nearest p points in N for each point in M
%
%p = N - 1;  %integer between m+1 and N
p = N-1;  %integer between m+1 and N

d_matrix = zeros(N,1);

for i=1:M

    for j=1:N

        d_matrix(j) = 0;

        for k=1:m

            d_matrix(j) = d_matrix(j) + norm(y_ref(i,k) - xi_ref(j,k));

        end

    end

    [z,index(i,:)] = sort(d_matrix);

    for j=1:p

        ip = (i-1)*p + j;

        points(ip,:) = xi(index(i,j),:);

    end

end

%
% formulate the least square problem
%
for np = 1 : M

    A = [1 points((np-1)*p+1,:)];

    for i = (np-1)*p+2 : np*p

        A = [A; 1 points(i,:)];

    end

    B = f(index(np,1));

    for i=2:p

        B = [B; f(index(np,i))];
    end

    z = A \ B;

    if np == 1

       b_matrix = z(2:m+1);

    else

       b_matrix = [b_matrix z(2:m+1)];

    end

end

%construct the covariance matrix

C = 0;

for i=1 : M

    z = b_matrix(:,i);

    C = C + z * z';

end

C=C/M;

[W D] = eig(C);

[lambda, idx] = sort(diag(D), 'descend');

W = W(:,idx);
eta1 = W(:,1);
eta2 = W(:,2);

% Eigenvalue Plot
figure(1);
semilogy(1:length(lambda),abs(lambda)./lambda(1),'ko','linewidth',2,'MarkerFaceColor','k');
xlabel('$$\mathrm{Index~(i)}$$','interpreter','latex','fontsize',20);
ylabel('$$\mathrm{Eigenvalue~(\lambda_i)}$$','interpreter','latex','fontsize',20);
set(gca,'TickLabelInterpreter','Latex','fontsize', 18);
%title('comparing dominant eigenvalues');
box on;
print -depsc eig_Zf3.eps

% univariate
%
figure(2);
g1 = eta1'*xi_ref';
plot(g1,f, 'ko', 'markerfacecolor', 'k')
set(gca, 'fontsize', 20);
xlabel('<eta1, x>');
ylabel('f(x)');
print -dpng ssp_Zf3.png

end

















