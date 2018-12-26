close all
clear all
rng(123);

n = 10;
x = linspace(-1,1,n);
y = linspace(-1,1,n);
[XX,YY] = meshgrid(x,y);
%p = 3.*XX+4.*YY;
p = XX.^2+YY.^2;
for j = 1:n
    for i = 1:n
        p(i,j) = p(i,j) + randn(1);
    end
end

Z = exp(-(p.^2)); % nxn


% PCA computation
sigma = Z(:)*Z(:)'; % n2xn2: n2 = n^2
[U,S,V] = svd(sigma); 
Ureduce = U(:,1); % n2x1
ze = Ureduce'*Z(:); % (1xn2)*(n2x1) = 1x1
Z_rec = ze*Ureduce'; % (1x1)*(1xn2) = 1xn2

figure(1);
pcolor(Z);
shading interp;
print -dpng origZ.png

Z_rec = reshape(Z_rec,n,n);
figure(2);
pcolor(Z_rec);
shading interp;
print -dpng recZ.png





