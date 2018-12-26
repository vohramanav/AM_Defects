ffffclose all;
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

nom = [vel,pow,pre,S,E,rho,C0_Cp,C1_Cp,C2_Cp,C0_k,C1_k,C2_k];

L = zeros(1,m); U = zeros(1,m);
L(1,:) = 0.9.*nom(1,:); U(1,:) = 1.1.*nom(1,:);

N = 20; % 20 samples in the 12D parameter space

% generate 20 samples in [-1,1]
xi  = -1 + 2*lhsdesign(N,m);
pts_x = gen_samples(xi,L,U);

function pts_x = gen_samples(xi,L,U)
pts_x = zeros(size(xi,1),size(xi,2));

for i = 1:size(pts_x,1)
  for j = 1:size(pts_x,2)
    pts_x(i,j) = L(1,j) + 0.5*(U(1,j)-L(1,j)).*(xi(i,j)+1);
  end
end

% Save physical points to a file
save('pts_gradfree_N20.txt','pts_x','-ASCII');
csvwrite('pts_N20.txt',pts_x)

end 
