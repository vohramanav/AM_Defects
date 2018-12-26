% visualizing the 2D stress data
close all
clear all

d = load('Stress_midSection_Dec20.mat');
x = zeros(1,32); y = zeros(1,14);
for i = 1:32
  x(1,32-i+1) = d.Stress_midSection((i-1)*10+1,1);
end
for i = 1:10
  y(1,10-i+1) = d.Stress_midSection(i,3);
end

for i = 1:4
y(1,14-i+1) = d.Stress_midSection(320+i,3);
end

[XX,YY] = meshgrid(x,y);

Z = zeros(14,32);

for i = 1:10
    for j = 1:32
        Z(10-i+1,32-j+1) = d.Stress_midSection((j-1)*10+i,5);
    end
end

for i = 1:4
    for j = 1:32
        Z(14-i+1,32-j+1) = d.Stress_midSection(320+(j-1)*4+i,5);
    end
end

save('Z_S11.txt','Z','-ASCII');

figure(1);
pcolor(Z);
shading interp;
colorbar();
print -dpng S11.png


