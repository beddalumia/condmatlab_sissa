function drawcrystal(x,y,z,a,r0,Lx,Ly,Lz,Nb)
N = length(x)/Nb;
% Displaying Crystal
figure("Name", "Crystal")
[Sx,Sy,Sz] = sphere;
for i = 1:Nb*N
    if mod(i,Nb) == 0
    if ~((x(i) >= Lx-a/2 || y(i) >= Ly-a/2 || z(i) >= Lz-a/2) && N <= 8)
    surf(x(i)+Sx*r0,y(i)+Sy*r0,z(i)+Sz*r0,'FaceColor',[0 .75 .75],'EdgeColor','none','FaceLighting','gouraud'); hold on
    end
    elseif mod(i,Nb) == 1
    if ~((x(i) >= Lx-a/2 || y(i) >= Ly-a/2 || z(i) >= Lz-a/2) && N <= 8)
        surf(x(i)+Sx*r0,y(i)+Sy*r0,z(i)+Sz*r0,'FaceColor',[.75 .75 0],'EdgeColor','none','FaceLighting','gouraud'); hold on
    end
    elseif mod(i,Nb) == 2
    if ~((x(i) >= Lx-a/2 || y(i) >= Ly-a/2 || z(i) >= Lz-a/2) && N <= 8)
        surf(x(i)+Sx*r0,y(i)+Sy*r0,z(i)+Sz*r0,'FaceColor',[.75 0 .75],'EdgeColor','none','FaceLighting','gouraud'); hold on
    end
    elseif mod(i,Nb) == 3
    if ~((x(i) >= Lx-a/2 || y(i) >= Ly-a/2 || z(i) >= Lz-a/2) && N <= 8)
        surf(x(i)+Sx*r0,y(i)+Sy*r0,z(i)+Sz*r0,'FaceColor',[.75 .75 .75],'EdgeColor','none','FaceLighting','gouraud'); hold on
    end
    end
end
plotcube([Lx Ly Lz],[-a/2 -a/2 -a/2],.1,[0.3 0.3 0.3]); % Conventional Cell Displying
axis equal
view(24.2576,19.5898);
lightangle(24.2576,19.5898); % add some light
end