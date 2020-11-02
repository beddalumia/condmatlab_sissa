a = 2.8; % Lattice parameter
ax = 3*a;
ay = sqrt(3)*a;
b_1 = [1,0]*a;               %%%%%%%%%%%%%%%%%%%%%%
b_2 = [-0.5, sqrt(3)/2]*a;   % Bond displacements %
b_3 = [-0.5, - sqrt(3)/2]*a; %%%%%%%%%%%%%%%%%%%%%%
eps = 0; %-10.7 eV to vacuum but since eps = EF and we want to set EF = 0...
t = 2.5; %eV So we'll keep the energies in electronvolt

M = 250;
ky_step = 2*pi/(ay*M);
ky_min= -pi/ay;
ky_max= +pi/ay;

Nmax=15;


for N = 2:1:Nmax
    
    k_y = [ky_min:ky_step:ky_max];
    eigs = zeros(length(k_y),4*N);
    
        for j = 1:length(k_y)
            ky = k_y(j);
            
            H_cell = zeros(4,4);
            H_cell(1,1) = eps;
            H_cell(2,2) = eps;
            H_cell(3,3) = eps;
            H_cell(4,4) = eps;
            H_cell(1,2) = -t*exp(1i*ky*b_1(2));
            H_cell(2,3) = -t*(exp(-1i*ky*b_2(2))+exp(-1i*ky*b_3(2)));
            H_cell(3,4) = -t*(exp(1i*ky*b_1(2)));
            H_cell(3,2) = conj(H_cell(2,3));
            H_cell(2,1) = conj(H_cell(1,2));
            H_cell(4,3) = conj(H_cell(3,4));
            
            H = zeros(N*4,N*4);

            for ii= 1:N
                H((ii-1)*4+1:(ii)*4, (ii-1)*4+1:(ii)*4) = H_cell;
            end

            PBC_Element = -t*(exp(1i*ky*b_2(2))+exp(1i*ky*b_3(2)));
            H(1,N*4) = PBC_Element;
            H(N*4,1) = conj(PBC_Element);

            for ii = 2:N
                H(4*(ii-1),4*(ii-1)+1)= PBC_Element;
                H(4*(ii-1)+1,4*(ii-1)) = conj(PBC_Element);

            end
            
            eigs(j,:) = eig(H);
    
        end
        
        figure(N) 
        %CM = lines(length(eigs(1,:)));
        for n = 1:length(eigs(1,:))
            plot(k_y,eigs(:,n),'black','LineWidth',1.5)%,'.-','color',CM(n,:))
            xlabel('k') 
            ylabel('E(k)') 
            xlim([ky_min,ky_max]);
            hold on
        end
        set(gcf,'papersize',[6 5])
        print(gcf,sprintf('Armchair1DBands%d',N),'-fillpage','-dpdf','-painters')
end