


data_N = load('number_particles.txt');
plot(data_N(:,1),data_N(:,2),'-','LineWidth',1.8)
hold on
plot(data_N(:,1),data_N(:,3),'-','LineWidth',1.8)
hold on
plot(data_N(:,1),data_N(:,4),'-','LineWidth',1.8)
xlabel('lattice')
legend('up','dn','total')
title('Number of particles')
grid on

%data_N = load('number_particles.txt');

data_corr = load('correlation.txt');
plot(data_corr(:,1),abs(data_corr(:,2)),'-','LineWidth',1.8)
hold on
plot(data_corr(:,1),abs(data_corr(:,3)),'-','LineWidth',1.8)

plot(data_corr(:,1),(data_corr(:,2)),'-','LineWidth',1.8)

grid on
semilogy(abs(data_corr(:,end)))
semilogy((data_corr(:,end)))

xlabel('lattice')
title('correlation function')
grid on

lattice = data_corr(:,1);
corr = log(abs(data_corr(:,2)));
%corr = log(data_corr(:,2));

data_corr = load('correlation_3.2.txt');
S = sparse(data_corr(:,1),data_corr(:,2),data_corr(:,3));
S = full(S);
Cmatrix = S + transpose(S) - diag(diag(S));
[V,D] = eig(Cmatrix);
plot(diag(D),'-','LineWidth',1.8)

plot(Cmatrix(3,3:end))

%Cmatrix = load('Cmatrix.txt');
%[V,D] = eig(Cmatrix);
V = V';
cmatrix_file = 'cmatix.mat';
save(cmatrix_file,'Cmatrix','D','V');

V1 = load('eigenvectors.txt');
V1 = V1';
V2 = load('eigenvectors_t0.txt');
for i = 1:20
    subplot(4,5,i);
    plot(V(i,:))
    %plot(V1(i,:))
    %plot(V2(i,:))
    a = num2str(D(i,i));
    if i == 1
        title(['eigenvalue: ',a])
    else
        title(a)
    end
    grid on
end

plot(diag(D),'-','LineWidth',1.8)

Hamiltonian1 = load('Hamiltonian1.txt');
Hamiltonian2 = load('Hamiltonian2.txt');

H1 = load('Hamiltonian1.txt');
eig(H1)
H2 = load('Hamiltonian2.txt');
eig(H2)
H12 = H1+H2;
eig(H12)

Hamiltonian_surf = load('Hamiltonian.txt');
surf(Hamiltonian_surf)

%%==========================extrapolation==================

extra1 = zeros(11,2);
for i = 1:11
    extra1(i,2) = (unnamed(i+1,3)-unnamed(1,3))/(i*10);
    extra1(i,1) = 1/((i+2)*10);
    %extra1(i,1) = (i+2)*10;
end
plot(extra1(:,1),extra1(:,2),'-','LineWidth',1.8)
xlabel('length: i*10-20')
ylabel('energy per site')
title('R=2.4,t1 = 0.1707,t2 = -0.03225,U = 0.84559659')
grid on
hold on
extra2 = zeros(5,2);
%extra2 = zeros(8,2);
for i = 1:5
    %unnamed(i,2) = (unnamed(1,i+7)-unnamed(1,i+3))/(40);
    %unnamed(i,1) = (i+1)*10; 
    %extra2(i,2) = (unnamed(i*2+4,3)-unnamed(i+1,3))/(30+i*10);
    extra2(i,2) = (unnamed(i*2+4)-unnamed(i+1))/(30+i*10);
    %extra2(i,2) = (unnamed(i*1+2)-unnamed(1))/(i*20+20);
    %extra2(i,1) = (i+1+4)*10; 
    extra2(i,1) = 1/((i+3)*10);
end
plot(extra2(:,1),extra2(:,2),'-','LineWidth',1.8)
tlattice = extra2(:,1);
tenergy = extra2(:,2);

xlabel('length: i*10 (50)')
ylabel('energy per site')
title('R=2.4,t1 = 0.1707,t2 = -0.03225,U = 0.84559659')
grid on

plot(extra1(4:11,1),extra1(4:11,2),'-','LineWidth',1.8)
hold on
plot(extra2(4:11,1),extra2(4:11,2),'-','LineWidth',1.8)
legend('method1: i*10-20','method2: i*10 (50)')
%%==========================
tlattice1 = zeros(5,1);
gap1 = zeros(5,1);
tlattice2 = zeros(5,1);
for i = 1:5
    tlattice2(i) = 1/tlattice1(i);
end
plot(tlattice2,gap1,'-','LineWidth',1.8)
%================

corr = zeros(1,20);
for n = 1:20
corre = 0;
    for i = 1:14
        corre = corre + 2*cos(i*((2*pi)/60)*n);
    end
corre = corre + 1 + cos((15)*((2*pi)/60)*n);
corre = corre*(1/60);
corr(n) = corre;
end


for i = 1:14
    corre = corre + 2*cos(i*((2*pi)/60)*1);
end
corre = corre + 1 + cos((15)*((2*pi)/60)*1);
corre = corre*(1/60);

%%==========================

data1 = load('Data_01_E-8.txt');
data2 = load('Data_015_E-8.txt');
data3 = load('Data_02_E-8.txt');

% plot(data1(:,1),data1(:,2),'-','LineWidth',1.8)
% xlabel('real time step')
% %legend('<Sz1>')
% hold on
% plot(data1(:,1),data1(:,3),'-','LineWidth',1.8)
% xlabel('real time step')
% legend('<Sz1>','<Sz2>')
% hold on
% grid on
plot(data1(:,1),data1(:,4),'-','LineWidth',1.8)
hold on
plot(data2(:,1),data2(:,4),'-','LineWidth',1.8)
hold on
plot(data3(:,1),data3(:,4),'-','LineWidth',1.8)
%xlabel('real time step')
%legend('SvN')
%plot(data1(:,1),data1(:,5),'-','LineWidth',1.8)
xlabel('real time step')
legend('0.01','0.15','0.2')
title('SvN of various tau')
grid on

plot(data1(:,1),data1(:,6),'-','LineWidth',1.8)
hold on
plot(data2(:,1),data2(:,6),'-','LineWidth',1.8)
hold on
plot(data3(:,1),data3(:,6),'-','LineWidth',1.8)
%xlabel('real time step')
%legend('SvN')
%plot(data1(:,1),data1(:,5),'-','LineWidth',1.8)
xlabel('real time step')
legend('0.01','0.15','0.2')
title('SvN of various tau')
grid on


plot(data1(:,1),data1(:,6),'-','LineWidth',1.8)
title('energy bond SS tau = 0.15')
xlabel('real time step')
grid on

fig = figure('Color','w');
%plotyy
%[ax, s1h1 s1h3] = plotyy(data1(:,1),data1(:,2:3),data1(:,1),data1(:,4),'plot');
[ax, s1h1 s1h3] = plotyy(data1(:,1),data1(:,4),data1(:,1),data1(:,5),'plot');
%sig1 color
sig1col = [0 150 150]/255;
%sig1log color
sig1logcol = [210 30 50]/255;
%style the plot
set(s1h1,'Color',sig1col,'LineWidth',3);
set(s1h3,'Color',sig1logcol,'LineWidth',3);
set(ax(1),'YColor',sig1col);
set(ax(2),'YColor',sig1logcol);

%x-axis and y-axis labels
xlabel('$x$','Interpreter','latex');
set(get(ax(1),'Ylabel'),'String','$\mbox{SvN}$','Interpreter','latex')
set(get(ax(2),'Ylabel'),'String','$\mbox{m}$','Interpreter','latex')
%add legend
%leg = legend('$\mbox{$Sz$}$ ','$\mbox{$Sz$}$ ','$\mbox{$SvN$}$ ', 'Location', 'NorthWest');
leg = legend('$\mbox{$SvN$}$ ','$\mbox{$m$}$ ', 'Location', 'NorthWest');

set(leg,'Interpreter','latex');

grid on
