a = 2.8; % Lattice parameter
ax = 3*a;
ay = sqrt(3)*a;
b_1 = [1,0]*a;               %%%%%%%%%%%%%%%%%%%%%%
b_2 = [-0.5, sqrt(3)/2]*a;   % Bond displacements %
b_3 = [-0.5, - sqrt(3)/2]*a; %%%%%%%%%%%%%%%%%%%%%%
eps = 0; %-10.7 eV to vacuum but since eps = EF and we want to set EF = 0...
t = 2.5; %eV So we'll keep the energies in electronvolt

M = 250;
kx_step = 2*pi/(ax*M);
kx_min= -pi/ax;
kx_max= +pi/ax;

Nmax=15;

for N = 2:1:Nmax
    
    figure(N)
    
    ky_step = 2*pi/(ay*N);
    if mod(N,2) == 0
    [k_y,k_x] = meshgrid(-N/2*ky_step:1*ky_step:N/2*ky_step, kx_min:kx_step:kx_max);
    else
    [k_y,k_x] = meshgrid(-(N-1)/2*ky_step:1*ky_step:(N-1)/2*ky_step, kx_min:kx_step:kx_max);    
    end
    eigs = zeros(length(k_x(:,1)),4);
    
    for i = 1:length(k_y(1,:))

        for j = 1:length(k_x(:,1))
            ky=k_y(1,i);
            kx=k_x(j,1);
                        
            H_cell = zeros(4,4);
            H_cell(1,1) = eps;
            H_cell(2,2) = eps;
            H_cell(3,3) = eps;
            H_cell(4,4) = eps;
            H_cell(1,2) = -t*exp(-1i*(kx*b_3(1)+ky*b_3(2)));
            H_cell(1,4) = -t*exp(-1i*(kx*b_1(1)+ky*b_1(2)));
            H_cell(2,3) = -t*exp(+1i*(kx*b_1(1)+ky*b_1(2)));
            H_cell(3,4) = -t*exp(-1i*(kx*b_2(1)+ky*b_2(2)));
            H_cell(3,2) = conj(H_cell(2,3));
            H_cell(2,1) = conj(H_cell(1,2));
            H_cell(4,3) = conj(H_cell(3,4));
            H_cell(4,1) = conj(H_cell(1,4));

            PBC_Element25 = -t*exp(+1i*(kx*b_2(1)+ky*b_2(2)));
            PBC_Element38 = -t*exp(-1i*(kx*b_3(1)+ky*b_3(2)));

            H_cell(1,2) = H_cell(1,2) + conj(PBC_Element25);
            H_cell(2,1) = conj(H_cell(1,2));
            
            H_cell(3,4) = H_cell(3,4) + PBC_Element38;
            H_cell(4,3) = conj(H_cell(3,4));
                                   
            eigs(j,:) = eig(H_cell);

        end
        
        CM = lines(length(k_y(1,:)));
        for jj = 1:length(eigs(j,:))
            plot(k_x,eigs(:,jj),'-','color',CM(i,:),'LineWidth',1.5)
            xlabel('k') 
            ylabel('E(k)')
            xlim([kx_min,kx_max]);
            hold on
        end

    end
   
    set(gcf,'papersize',[6 5])
    print(gcf,sprintf('ZigZag2DBands%d',N),'-fillpage','-dpdf','-painters')
   
end