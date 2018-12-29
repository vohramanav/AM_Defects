% visualizing the 2D stress data
close all
clear all

d = load('../Data/allStressDec25.mat');
x = zeros(1,32); y = zeros(1,14);
S = d.allStress(:,:,1);

for i = 1:32
  x(1,32-i+1) = S((i-1)*10+1,1);
end
for i = 1:10
  y(1,10-i+1) = S(i,3);
end

for i = 1:4
y(1,14-i+1) = S(320+i,3);
end

[XX,YY] = meshgrid(x,y);

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

save('Z_S11_sam1.txt','Z','-ASCII');

%figure(1);
%pcolor(Z);
%shading interp;
%axis equal;
%xlim([1,32]);
%ylim([1,14]);
%colorbar();
%print -dpng S11_sam15.png
