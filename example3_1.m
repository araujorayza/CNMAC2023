clear all;
clc;
close all;

%system
A{1}=[-10 30;0 -20];
A{2}=[-10 -30;0 -20];
n=2;

h{1} = @(x1,x2) (x1.^2+x2.^2)/50;
h{2} = @(x1,x2) 1-(x1.^2+x2.^2)/50;

dh{1} = @(x1,x2) [ 2.*x1/50,  2.*x2/50];
dh{2} = @(x1,x2) [-2.*x1/50, -2.*x2/50];

G=[1];
b=0.35;
l=0.2;

Rset=1:n;

Z=-5:0.01:5;

%LMI calculations
LMIS=[];
for j=G
    P{j} = sdpvar(n,n,'symmetric');
    R{j} = sdpvar(n,n,'full');
    L{j} = sdpvar(n,n,'full');
end

for j=Rset
    for k=G
        Upsilon{k,j} = [L{k}*A{j}+A{j}'*L{k}',   (P{k}-L{k}'+R{k}*A{j})';
                        P{k}-L{k}'+R{k}*A{j},        -R{k}-R{k}'];
        LMIS = [LMIS, Upsilon{k,j} <= 0];
    end
end


opts=sdpsettings;
opts.solver='sedumi';
opts.verbose=0;

sol = solvesdp(LMIS,[],opts);
p=min(checkset(LMIS));
if p > 0
    for k = G
        P{k} = double(P{k})
        R{k} = double(R{k})
        L{k} = double(L{k})
    end
else
    display('Infeasible')
    P=[];
end

% Set estimation
V = @(x1,x2) sum(arrayfun(@(k) [x1;x2]'*h{k}(x1,x2)*P{k}*[x1;x2],G));
hdot = @(x1,x2,k) sum(arrayfun(@(j) dh{k}(x1,x2)*h{j}(x1,x2)*A{j}*[x1;x2],Rset));
Dset = @(x1,x2) sum(arrayfun(@(k) [x1;x2]'*hdot(x1,x2,k)*P{k}*[x1;x2],G));

figure; %this is where you obtain the information about b
fs=fsurf(V);
title('V(x) along the system trajectories inside Z')

figure; %this is where you obtain the information about l
fc=fcontour(V);
fc.LevelList = linspace(0,b,10);
fc.DisplayName = strcat('V level sets: ',num2str(fc.LevelList));
hold on;
fc=fcontour(Dset);
fc.LineColor='k';
fc.DisplayName='D';
fc.LevelList=[0];
title('V level-sets and D');
legend;
colorbar;
caxis([0,b]);

figure;
fc=fcontour(V);
fc.LineColor ='r';
fc.LineStyle ='--';
fc.LineWidth = 2;
fc.DisplayName=num2str(l);
fc.LevelList=[l];

hold on

fc=fcontour(V);
fc.LineColor ='b';
fc.LineStyle =':';
fc.LineWidth = 2;
fc.DisplayName=num2str(b);
fc.LevelList=[b];

fc=fcontour(Dset);
fc.LineColor ='k';
fc.LineStyle ='-.';
fc.LineWidth = 2;
fc.DisplayName = 'SetD';
fc.LevelList=[0];

hold off
title('V level-sets and D');