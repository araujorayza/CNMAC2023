%% This sections creates the drone model and calculates the LMI results
clear all; close all; clc;
ControlType = 'CNMAC2023';
Modeltype = 1;
gamma = [3.8195    0.1614    3.8117    0.3855    1.6844    1.7983    4.4994    2.0476];
global small_k 
small_k = 0.1; %gain value of the error dynamic gain
Saturation = false;
ControlDesign; 
% available variables at this point:
% A = system matrix (drone+controller)
% B, C = model matrices (drone+controller)
% G = conjunto G do artigo CNMAC
% K = ganho do controlador calculado (por uma LMI do mozelli que ja usei
% antes)
% P = matrizes da funçao V calculadas na LMI
% R = matrix auxiliar
% L = matrix auxiliar
% h = funcoes de pertinencia do modelo

%% This section was supposed to calculate the sets
% Set estimation
V = @(x1,x2) sum(arrayfun(@(k) [x1;x2]'*h{k}(x1,x2)*P{k}*[x1;x2],G));
hdot = @(x1,x2,k) sum(arrayfun(@(j) dh{k}(x1,x2)*h{j}(x1,x2)*A{j}*[x1;x2],Rset));
Dset = @(x1,x2) sum(arrayfun(@(k) [x1;x2]'*hdot(x1,x2,k)*P{k}*[x1;x2],G));

%calculate V and D
meshPoints=500;
tol=10/meshPoints;
x = linspace(-5,5,meshPoints);
y = linspace(-5,5,meshPoints);
[X,Y]=meshgrid(x,y);
for i=1:length(x)
    for j = 1:length(y)
        Ve(i,j) = V(X(i,j),Y(i,j));
        De(i,j) = Dset(X(i,j),Y(i,j));
    end
end

%calculate b
b=min([min(Ve(:,1)), min(Ve(:,end)), min(Ve(1,:)), min(Ve(end,:))])
b=fix(b*1e2)/1e2;

figure(1);
[~,c]=contour(X,Y,Ve,linspace(0,b,5),'r','ShowText','on','DisplayName','V')
hold on
[~,d]=contour(X,Y,De,[0,fix(max(max(De))*1e2)/1e2],'b','ShowText','on','DisplayName','D')
legend;
%from the graph
l = 0.18;


