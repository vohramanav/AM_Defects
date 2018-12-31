close all
clear all
rng(123);

Stot = 20;
K = 20; % choose only 20 components to reconstruct the data
%extract_stress(Stot);
[X,Ureduce] = compute_Ureduce(Stot,K);
recover(Stot,X,Ureduce,K);

function data = extract_stress(Stot)
% Extract stress fields for the 20 sample points
d = load('../Data/allStressDec25.mat');

for nsamp = 1:Stot
  S = d.allStress(:,:,nsamp);
  Z = zeros(14,32);

  for i = 1:10
      for j = 1:32
          Z(10-i+1,32-j+1) = S((j-1)*10+i,5);
      end
  end
  
  for i = 1:4
      for j = 1:32
          Z(14-i+1,32-j+1) = S(320+(j-1)*4+i,5);
      end
  end

  fname = sprintf('Z_files/Z_S11_sam%d.txt',nsamp);
  save(fname,'Z','-ASCII');
end

end

function [X,Ureduce] = compute_Ureduce(Stot,K)
r = 14; c = 32; % number of rows and columns in Z
m = Stot; n = r*c;
X = zeros(m,n);

% Computing the overall dominant eigenvector
for i = 1:Stot
  fname = sprintf('Z_files/Z_S11_sam%d.txt',i);
  Z = load(fname);
  X(i,:) = Z(:); % mxn: n = r*c
end
sigma = (1./Stot).*(X'*X); % nxn 
[U,S,V] = svd(sigma); 
Ureduce = U(:,1:K); % nxK, the overall dominant eigenvector
save('Ureduce.txt','Ureduce','-ASCII');

end

function r = recover(Stot,X,Ureduce,K)
r = 14; c = 32; % number of rows and columns in Z
%X_rec = Zproj*Ureduce'; % (mxK)*(K*n) = (mxn)

max_diff_K = zeros(K,1); % store max val of max_diff_Ki for each K
for Ki = 1:K
  Z_proj = X*Ureduce(:,1:Ki); % (mxn)*(n*K) = (mxK)
  X_rec = Z_proj*Ureduce(:,1:Ki)';
  max_diff_Ki = zeros(Stot,1); % store max difference for a ith K i.e Ki
  for j = 1:Stot
      max_diff_Ki(j,1) = max(X(j,:)-X_rec(j,:));
  end
  max_diff_K(Ki,1) = max(max_diff_Ki);
end

figure(1)
plot(1:K,max_diff_K,'--bo','MarkerFaceColor','b','MarkerSize',3);
xlabel('Number of Principal Components (K)','Fontsize',20);
ylabel('Error','Fontsize',20);
set(gca,'fontsize',18);
print -dpng error_k.png

% K=3 seems like a good choice
ko = 3; % ko - K optimal
Zf = X*Ureduce(:,1:ko);
save('Zf.txt','Zf','-ASCII');

Z_proj = load('SZF.txt');
X_rec = Z_proj*Ureduce(:,1:3)';

X1 = reshape(X(20,:),r,c);
figure(1);
pcolor(X1);
shading interp;
axis equal;
xlim([1,32]);
ylim([1,14]);
colorbar();
print -dpng origZ_sam20.png

X_rec1 = reshape(X_rec(1,:),r,c);
figure(2);
pcolor(X_rec1);
shading interp;
axis equal;
xlim([1,32]);
ylim([1,14]);
colorbar();
print -dpng recZ_sam20.png

end
