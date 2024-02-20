%Script to calculate the controlller gains for the bebop quadrotor
%mostly based on various T-S Fuzzy models
[A,B,h,C] = ErrorModeling(Modeltype,gamma)

switch ControlType
    case 'openloop'
        K = zeros(4);
        disp('Virtual controller set to open loop');

    case 'MozelliTeo6'
        fi=0.01*ones(size(A,2));
        mu=0.1;

        if(Modeltype == 2)
            K = LMI_Teo6MozelliMOD(A,B,fi,mu);
        else
            K = LMI_Teo6Mozelli(A,B,fi,mu);
        end
        ControlType = 'PDC';
    case 'SereniTeo2'
        %for i=0:1:100
        K = SereniTeo2(A,B,C,0.06);
        %end
        ControlType = 'PDC_SOF';
    case 'SereniTeo2_Ulim'
        K = SereniTeo2_Ulim(A,B,C,0.04,2,0.8);
        ControlType = 'PDC_SOF';
    case 'WeiTeo1'
        r = 0.21; % Valor m�ximo para os autovalores de W
        K = WeiTeo1(A,B,C,r);
        ControlType = 'PDC_SOF';
    case 'teste'
        mu=0.01;
        K = lmi_Dong2013(A,B,C,mu)
        ControlType = 'PDC_SOF';
    case 'ICUAS'
        MaxRotSpd
        %         phi=0.95;
        %         Mu=1;
        %         K = ICUAS_Teo1(Modeltype,A,B,phi,Mu);
        phi = 0.95; Mu = 1; ga = 2;
        K = ICUAS_Teo2(Modeltype,A,B,B,C,phi,Mu,ga);
        ControlType = 'PDC';
    case 'PhDpaper'
        MaxRotSpd
        %         phi=0.95;
        %         Mu=1;
        %         K = ICUAS_Teo1(Modeltype,A,B,phi,Mu);
        phi = 0.95; Mu = 1; ga = 2;
        K = ICUAS_Teo2(Modeltype,A,B,B,C,phi,Mu,ga);
        ControlType = 'PDC';
    case 'CNMAC2023'
        [K,P,R,L,A,G,Rset] = CNMAC2023(Modeltype,A,B);
        ControlType = 'PDC';
    otherwise
        K=[];
        disp('The controller you chose is not an option!')
end

SimStruct.sys.sat = Saturation;

SimStruct.controller.type = ControlType;
SimStruct.controller.model.type = Modeltype;
SimStruct.controller.model.A = A;
SimStruct.controller.model.B = B;
SimStruct.controller.model.C = C;
SimStruct.controller.model.h = h;
SimStruct.controller.gain    = K;

%% auxiliar functions
function F = LMI_Teo6Mozelli(A,B,fi,mu)
% This implements Theorem 6 on
% "A systematic approach to improve multiple Lyapunov function stability
% and stabilization conditions for fuzzy systems" by
% L. A. Mozelli, R. M. Palhares, G. S. C. Avellar. 2009
% -------------------------------------------------------------------
r = size(B,2);
nA = size(A{1},2);
mB = size(B{1},2);

T = cell(1,r);
Y = sdpvar(nA,nA,'symmetric');
R = sdpvar(nA,nA,'full');
S = cell(1,r);

Tfi = 0;
for i = 1:r
    T{i} = sdpvar(nA,nA,'symmetric');
    S{i} = sdpvar(nA,mB);
    Tfi = Tfi + fi(i)*(T{i}+Y);
end

Csi = cell(r);
Csibarra = cell(r);
for i = 1:r
    for j = 1:r
        star = T{i} - mu*(A{i}*R'-B{i}*S{i}') + R;
        Csi{i,j} = [Tfi-A{i}*R'-R*A{i}'+B{i}*S{i}'+S{i}*B{i}'  star';
            star                            mu*(R+R')];
    end
end

for i = 1:r
    for j = 1:r
        Csibarra{i,j} = Csi{i,j} + Csi{j,i};
    end
end

Restr = [];
for j = 1:r
    for i = 1:j
        Restr = [Restr T{i}>=0 ...
            T{i}+Y>=0 ...
            Csi{i,i} <=0 ...
            Csibarra{i,j}<=0];
    end
end
% Configurando o Solver.
opts=sdpsettings;
opts.verbose=0;

% Resolvendo as LMIs
sol = solvesdp(Restr,[],opts);

% Verifica se as LMI's sao factiveis
p=min(checkset(Restr));
if (p > 0)
    R = double(R);
    Y = double(Y);
    F = cell(1,r);
    for i = 1:r
        T{i} = double(T{i});
        S{i} = double(S{i});
        F{i} = S{i}'/(R');
    end
else
    F=[];
    disp('LMIs Infactiveis')
end
end
%%
function K = LMI_Teo6MozelliMOD(A,B,fi,mu)
% Esta fun��o encontra condi��es suficientes para o projeto de
% controladores fuzzy, considerando realimenta��o de estados.
% Esse trabalho implementa o resultado do Teorema 6 de:
% A systematic approach to improve multiple Lyapunov function stability
% and stabilization conditions for fuzzy systems
% L. A. Mozelli, R. M. Palhares, G. S. C. Avellar. 2009
% -------------------------------------------------------------------
ri = size(B,2);
nA = size(A{1,1},2);     % Determina a dimens�o da matriz A
mB = size(B{1,1},2);     % Determina o n�mero de colunas de B
% Cell que vai armazenar os controladores locais
T = cell(1,ri);
S = cell(1,ri);
% Criando as vari�veis matriciais
Y = sdpvar(nA,nA,'symmetric');  % Declaracao de Y nxn simetrica
R = sdpvar(nA,nA,'full');       % Declaracao de R nxn full
% Criando as vari�veis Ti, Si e Tfi
Tfi = 0;
Restr = [];
for t1 = 1:ri
    T{1,t1} = sdpvar(nA,nA,'symmetric');  % Declaracao de Ti nxn simetrica
    S{1,t1} = sdpvar(nA,mB);              % Declaracao de Si nxm
    Tfi = Tfi + fi(1,t1)*(T{1,t1}+Y);
end
% % Criando o conjunto de restri��es LMI
for t1 = 1:ri
    Restr = [Restr T{1,t1}>=0 T{1,t1}+Y>=0 ...
        [Tfi-A{t1,t1}*R'-R*A{t1,t1}'+B{1,t1}*S{1,t1}'+S{1,t1}*B{1,t1}'  ...
        (T{1,t1}-mu*(A{t1,t1}*R'-B{1,t1}*S{1,t1}')+R)';
        T{1,t1}-mu*(A{t1,t1}*R'-B{1,t1}*S{1,t1}')+R  mu*(R+R')]<=0];
    for t2 = t1+1:ri
        Restr = [Restr ...
            [Tfi-A{t1,t2}*R'-R*A{t1,t2}'+B{1,t1}*S{1,t2}'+S{1,t2}*B{1,t1}'  ...
            (T{1,t1}-mu*(A{t1,t2}*R'-B{1,t1}*S{1,t2}')+R)';
            T{1,t1}-mu*(A{t1,t2}*R'-B{1,t1}*S{1,t2}')+R  mu*(R+R')]+...
            [Tfi-A{t2,t1}*R'-R*A{t2,t1}'+B{1,t2}*S{1,t1}'+S{1,t1}*B{1,t2}'  ...
            (T{1,t2}-mu*(A{t2,t1}*R'-B{1,t2}*S{1,t1}')+R)';
            T{1,t2}-mu*(A{t2,t1}*R'-B{1,t2}*S{1,t1}')+R  mu*(R+R')]<=0];
    end
end
% Configurando o Solver.
opts=sdpsettings;
% opts.solver='lmilab';
opts.verbose=0;
% Resolvendo as LMIs
sol = solvesdp(Restr,[],opts);
% Verifica se as LMI's sao factiveis
p=min(checkset(Restr));
if p > 0
    R = double(R);
    Y = double(Y);
    K = cell(1,ri);
    for t1 = 1:ri,
        T{1,t1} = double(T{1,t1});
        S{1,t1} = double(S{1,t1});
        K{1,t1} = S{1,t1}'*inv(R');
    end
else
    K=[];
    error('LMIs Infactiveis')
end
end
%%
function L = SereniTeo2(A,B,C,gamma)
N = size(A,2);
nA = size(A{1,1},2);
mB = size(B{1,1},2);
mC = size(C,1);

F = cell(1,N);
G = cell(1,N);
J = cell(1,N);
P = sdpvar(nA,nA,'symmetric');
H = sdpvar(mB,mB,'full');
for i=1:N
    F{i} = sdpvar(nA,nA,'full');
    G{i} = sdpvar(nA,nA,'full');
    J{i} = sdpvar(mB,mC,'full');
end
K = Quadratic(A,B,gamma);
LMIs = [P>=0];
for i=1:N
    a11 = A{i}'*F{i}' + F{i}*A{i} + K'*B{i}'*F{i}' + F{i}*B{i}*K + 2*gamma*P;
    a21 = P - F{i}' + G{i}*A{i} + G{i}*B{i}*K;
    a31 = B{i}'*F{i}' + J{i}*C - H*K;
    a22 = -G{i}-G{i}';
    a32 = B{i}'*G{i}';
    a33 = -H-H';

    LMIs = [LMIs,   [a11 a21' a31';
        a21 a22  a32';
        a31 a32  a33]<=0];
    for j = i+1:N
        a11 = A{i}'*F{j}' + K'*B{i}'*F{j}' + A{j}'*F{i}' + K'*B{j}'*F{i}';
        a11 = a11 + a11' + 4*gamma*P;

        a21 = 2*P - F{i}' -F{j}' + G{i}*A{j} + G{i}*B{j}*K + G{j}*A{i} + G{j}*B{i}*K;
        a31 = B{i}'*F{j}' + J{j}*C + J{i}*C + B{j}'*F{i}' - 2*H*K;
        a22 = -G{i}-G{i}'-G{j}-G{j}';
        a32 = B{i}'*G{j}' + B{j}'*G{i}';
        a33 = -2*H-2*H';

        LMIs = [LMIs,[  a11 a21' a31';
            a21 a22  a32';
            a31 a32  a33]<=0];
    end
end
opts=sdpsettings;
sol = optimize(LMIs,[],opts);
che=min(check(LMIs));
if che > 0
    disp('funcionou!')
    H = value(H);
    for i=1:N
        J{i} = value(J{i});
        P  = value(P);
        L{i} = -inv(H)*J{i};
    end

else
    L=[];
    error('nao funcionou')
end
end
%%
function K = Quadratic(A,B,beta)
N = size(A,2);
nA = size(A{1,1},2);
mB = size(B{1,1},2);
W = sdpvar(nA,nA,'symmetric');
Z = sdpvar(mB,nA,'full');

LMIs = [W>=0];

for i=1:N
    LMIs = [LMIs, A{i}*W + W*A{i}' + B{i}*Z + Z'*B{i}' + 2*beta*W <= 0];
end
opts=sdpsettings;
% opts.solver = 'lmilab';
sol = optimize(LMIs,[],opts);
che=min(check(LMIs));
if che > 0
    disp('quadratico funcionou')
    W = value(W);
    K = value(Z)*inv(W)
else
    K = []
    error('quadratico nao funcionou')
end
end
%%
function L = SereniTeo2_Ulim(A,B,C,gamma,ro,xb)
N = size(A,2);
nA = size(A{1,1},2);
mB = size(B{1,1},2);
mC = size(C,1);

F = cell(1,N);
G = cell(1,N);
J = cell(1,N);
K = cell(1,N);
P = cell(1,N);
H = sdpvar(mB,mB,'full');
for i1=1:N
    F{i1} = sdpvar(nA,nA,'full');
    G{i1} = sdpvar(nA,nA,'full');
    J{i1} = sdpvar(mB,mC,'full');
    P{i1} = sdpvar(nA,nA,'symmetric');
end
K = Quadratic_PDC(A,B,gamma);
T = cell(N,N,N);
for i1=1:N
    for i2 = 1:N
        for i3 = 1:N
            a11 = A{i1}'*F{i2}' + F{i2}*A{i1} + K{i3}'*B{i1}'*F{i2}' + F{i2}*B{i1}*K{i3} + 2*gamma*P{i2};
            a21 = P{i2} - F{i2}' + G{i2}*A{i1} + G{i2}*B{i1}*K{i3};
            a31 = B{i1}'*F{i2}' + J{i1}*C - H*K{i3};
            a22 = -G{i2}-G{i2}';
            a32 = B{i1}'*G{i2}';
            a33 = -H-H';
            T{i1,i2,i3}=[a11 a21' a31';
                a21 a22  a32';
                a31 a32  a33];
        end
    end
end
clear LMIs
LMIs = [-ro^2*eye(mB) H';H -eye(mB)]<=0;
for i1=1:N
    LMIs = [LMIs P{i1}>=0 T{i1,i1,i1}<=0];
    % ******************************
    % LMI que limita o u
    a11=-(ro/xb)^2*eye(nA);
    a21=J{i1}*C;
    a22=-eye(mB);
    LMIs = [LMIs, [a11 a21';a21 a22]<=0];
    % ******************************
    for i2 = 1:N
        if i1~=i2
            LMIs = [LMIs T{i1,i1,i2}+T{i1,i2,i1}+T{i2,i1,i1}<=0];
        elseif i1<i2
            for i3 = i2+1:N
                LMIs = [LMIs T{i1,i2,i3}+T{i3,i2,i1}+...
                    T{i2,i1,i3}+T{i3,i1,i2}+...
                    T{i1,i3,i2}+T{i2,i3,i1}<=0];
            end
        end
    end
end

opts=sdpsettings;
% opts.solver = 'lmilab'; % Rodei usando o Mozek
opts.verbose = 0;
sol = optimize(LMIs,[],opts);
che=min(check(LMIs));
if che > 0
    disp('funcionou!')
    H = value(H);
    for i1=1:N
        J{i1} = value(J{i1});
        P{i1} = value(P{i1});
        L{i1} = -H\J{i1};
    end

else
    L=[];
    disp('nao funcionou')
end
end
%%
function K = Quadratic_PDC(A,B,gamma)
% beta = 1;
N = size(A,2);
nA = size(A{1,1},2);
mB = size(B{1,1},2);
W = sdpvar(nA,nA,'symmetric');
Z = cell(1,N);
for i1=1:N
    Z{i1} = sdpvar(mB,nA,'full');
end
clear LMIs
LMIs = W>=0;
for i1=1:N
    tii=A{i1}*W+W*A{i1}'+B{i1}*Z{i1}+Z{i1}'*B{i1}'+2*gamma*W;
    LMIs = [LMIs tii <= 0];
    for i2=1:N
        if i1~=i2
            LMIs = [LMIs 2*tii/(N-1)+...
                A{i1}*W+W*A{i1}'+B{i1}*Z{i2}+Z{i2}'*B{i1}'+2*gamma*W+...
                A{i2}*W+W*A{i2}'+B{i2}*Z{i1}+Z{i1}'*B{i2}'+2*gamma*W<=0];
        end
    end
end
opts=sdpsettings;
% opts.solver = 'lmilab';
opts.verbose = 0;
sol = optimize(LMIs,[],opts);
che=min(check(LMIs));
if che > 0
    disp('PDC funcionou')
    W = value(W);
    for i1=1:N
        K{i1} = value(Z{i1})*inv(W);
    end
else
    K = [];
    error('PDC nao funcionou')
end
end
%%
function K = WeiTeo1(A,B,C,r)
N = size(A,2);
nA = size(A{1,1},2);
mB = size(B{1,1},2);
W = sdpvar(nA,nA,'symmetric');
Z = sdpvar(mB,nA,'full');
LMIs = W>=r*eye(nA);
for i1=1:N
    a11 = A{i1}*W+W*A{i1}'-B{i1}*Z-Z'*B{i1}';
    LMIs = [LMIs, [a11 Z';Z -eye(mB)]<=0];
end
LMIs = [LMIs, [-W/r Z';Z -eye(mB)]<=0];
% Configurando o Solver.
opts=sdpsettings;
% opts.solver='lmilab';
% opts.solver='sedumi';
opts.verbose=0;
% Resolvendo as LMIs
sol = optimize(LMIs,[],opts);
che=min(check(LMIs));
if che > 0
    display('WEI2018 - SIM')
    % Encontra o valor num�rico das matrizes
    W = value(W);
    P = inv(W);
    Z = value(Z);
    K = cell(1,N);
    %     AutoWZ=[]; AutoCW=[]; AutoA=[];
    for i1=1:N
        K{i1} = -Z*inv(W)*pinv(C);
        %         AutoWZ=[AutoWZ real(eig(A{i1}*W-B{i1}*Z))];
        %         AutoCW=[AutoCW real(eig((A{i1}-B{i1}*K{i1}*C)*W))];
        %         AutoA =[AutoA  real(eig(A{i1}-B{i1}*K{i1}*C))];
    end
    Z
    Y1=K{1}*C*W
    %     AutoWZ
    %     AutoCW
    %     AutoA
    %celldisp(K)
    %     MAX_Auto=max(max(real(AutoA)))
else
    K=[];
    error('infactivel')
end
end
%%
function L = lmi_Dong2013(A,B,C,mu)
rp = size(A,2);
% Condi��es de estabilidade baseadas no Dong 2013 com
nA = size(A{1,1},2);
nB = size(B{1,1},2);
nC = size(C,1);

% Criando as vari�veis matriciais
G = sdpvar(nC,nC,'full');
T = inv(C*C'); % Matriz definida em (9) de [J. Dong and G. Yang, 2013]
Q = cell(1,rp);
L = sdpvar(nB,nC,'full');
Restr = [];
for t1=1:rp
    Q{1,t1} = sdpvar(nA,nA,'symmetric');
    Restr = [Restr,Q{1,t1}>=0];
end
% Criando as estruturas matriciais usadas nas LMIs
LMI = cell(rp,rp);
for t1=1:rp
    for t2=1:rp
        a11=A{1,t1}*Q{1,t2}+Q{1,t2}*A{1,t1}'+...
            B{1,t1}*L*T*C+(B{1,t1}*L*T*C)';
        a21=C*Q{1,t1}-G*T*C+mu*(B{1,t1}*L)';
        a22=-mu*(G+G');
        LMI{t1,t2} = [a11 a21';a21 a22];
    end
end
% Criando as LMIs
for t1=1:rp
    for t2=t1:rp
        Restr = [Restr,LMI{t1,t2}+LMI{t2,t1}<=0];
    end
end
% Configurando o Solver.
opts=sdpsettings;
opts.savesolverinput=1;
opts.savesolveroutput=1;
% opts.solver='lmilab';
% opts.solver='sedumi';
opts.verbose=0;
% Resolvendo as LMIs
sol = optimize(Restr,[],opts);
che=min(check(Restr));
if che > 0
    disp('DONG2013 - SIM')
    % Encontra o valor num�rico das matrizes
    G = value(G);
    L = value(L)*inv(G);
    for i1=1:rp,
        Q{1,i1}  = value(Q{1,i1});
    end
else
    L=[];
    disp('DONG2013 - N�O')
end
end
%%
function K = ICUAS_Teo1(Modeltype,A,B,phi,Mu)
%LMIs implemeted based on Theorem 1 from paper "------" ICUAS 2020
theta1=0.35; theta2=-0.5; dec=0.25; % maxU = 0.8632    0.5923    7.3816    0.0005
[nx,nu]=size(B{1});
r = length(A);
W = sdpvar(nx,nx,'full');
Q = cell(1,r);
Y = cell(r);
LMIs = [];
for i=1:r
    Q{i}= sdpvar(nx,nx,'symmetric');
    Y{i}= sdpvar(nu,nx,'full');
    LMIs = [LMIs, Q{i}>=0];
end
if(Modeltype == 3 || Modeltype == 4)
    warning('Local Models not based on nonlinearity, check how the matrix G is made')
end
G=PdotConvex(log2(r));
eta = length(G);
Qbar = cell(1,eta);
for l=1:eta
    Qbar{l} = 0;
    for i=1:r
        Qbar{l} = Qbar{l} + phi*G(i,l)*Q{i};
    end
end

if(Modeltype == 2) %model with two indices
    B=B{1,1};
    for i = 1:r
        for j = 1:r
            for l = 1:eta
                a11 = Qbar{l} - A{i,j}*W - ...
                    W'*A{i,j}' + B*Y{i} + Y{i}'*B'+2*dec*Q{i};
                a21 = Q{i} + W' - Mu*(A{i,j}*W - B*Y{i});
                a22 = Mu*(W+W');
                Upsilon{i,j,l} = [a11, a21';
                    a21, a22];
            end
        end
    end
else
    B=B{1,1};
    for i = 1:r
        for l = 1:eta
            a11 = Qbar{l} - A{i}*W - ...
                W'*A{i}' + B*Y{i} + Y{i}'*B'+2*dec*Q{i};
            a21 = Q{i} + W' - Mu*(A{i}*W - B*Y{i});
            a22 = Mu*(W+W');
            Upsilon{i,l} = [a11, a21';
                a21, a22];
        end
    end
end


% Criando as LMIs
LMIs = [LMIs, W+W'<=theta2*eye(nx)];
for i=1:r
    SizeConstr = [theta1*eye(nu), Y{i};
        Y{i}',    eye(nx)];
    LMIs = [LMIs, SizeConstr>=0];
    for j=1:r
        if(i==j)
            continue;
        end
        for l=1:eta
            if(Modeltype == 2) %model with two indices
                LMIs = [LMIs, Upsilon{i,i,l}<=0];
                LMIs = [LMIs, 2/(r-1)*Upsilon{i,i,l} + ...
                    Upsilon{i,j,l} + Upsilon{j,i,l}<=0];
            else
                LMIs = [LMIs, Upsilon{i,l}<=0];
            end

        end
    end
end

% Configurando o Solver.
opts=sdpsettings;
% opts.solver='lmilab';
opts.solver='sedumi';
opts.verbose=0;
% Resolvendo as LMIs
sol = solvesdp(LMIs,[],opts);
che=min(checkset(LMIs));
if che > 0
    % Encontra o valor num�rico das matrizes
    W=double(W);
    K = cell(1,r);
    for i=1:r
        K{i} = double(Y{i})/W;
    end
else
    error('LMIs infact�veis')
    Q=[]; W=[]; K=[];
end
end
%%
function K = ICUAS_Teo2(Modeltype,A,B,Bw,C,phi,Mu,ga)
%LMIs implemeted based on Theorem 3 from paper http://dx.doi.org/10.1016/j.apm.2014.03.034
theta1=0.1; theta2=0.5; dec=0.01; % maxU = 0.8632    0.5923    7.3816    0.0005
theta1=0.5; theta2=0.5; dec=0.01;
[nx,nu]=size(B{1});
nc=size(C,1);
r = length(A);
W = sdpvar(nx,nx,'full');
Q = cell(1,r);
Y = cell(r);
LMIs = [];
for i=1:r
    Q{i}= sdpvar(nx,nx,'symmetric');
    Y{i}= sdpvar(nu,nx,'full');
    LMIs = [LMIs, Q{i}>=0];
end
if(Modeltype == 3 || Modeltype == 4)
    warning('Local Models not based on nonlinearity, check how the matrix G is made')
end
G=PdotConvex(log2(r));
eta = length(G);
Qbar = cell(1,eta);
for l=1:eta
    Qbar{l} = 0;
    for i=1:r
        Qbar{l} = Qbar{l} + phi*G(i,l)*Q{i};
    end
end

if(Modeltype == 2) %model with two indices
    B=B{1,1}; Bw=Bw{1,1};
    for i = 1:r
        for j = 1:r
            for l = 1:eta
                a11 = Qbar{l} + A{i,j}*W + ...
                    W'*A{i,j}' - B*Y{i} - Y{i}'*B'+2*dec*Q{i};
                a21 = Q{i} - W' + Mu*(A{i,j}*W - B*Y{i});
                a22 = -Mu*(W+W');
                a31=Bw';
                a32=Mu*a31;
                a33=-eye(nu)*(ga^2);
                a41=C*W;
                a42=zeros(nc,nx);
                a43=zeros(nc,nu);
                a44=-eye(nc);
                Upsilon{i,l}=[ a11 a21' a31' a41'
                    a21 a22  a32' a42'
                    a31 a32  a33  a43'
                    a41 a42  a43  a44];
            end
        end
    end
else
    B=B{1,1}; Bw=Bw{1,1};
    for i = 1:r
        for l = 1:eta
            a11 = Qbar{l} + A{i}*W + ...
                W'*A{i}' - B*Y{i} - Y{i}'*B'+2*dec*Q{i};
            a21 = Q{i} - W' + Mu*(A{i}*W - B*Y{i});
            a22 = -Mu*(W+W');
            a31=Bw';
            a32=Mu*a31;
            a33=-eye(nu)*(ga^2);
            a41=C*W;
            a42=zeros(nc,nx);
            a43=zeros(nc,nu);
            a44=-eye(nc);
            Upsilon{i,l}=[ a11 a21' a31' a41'
                a21 a22  a32' a42'
                a31 a32  a33  a43'
                a41 a42  a43  a44];
        end
    end
end


% Criando as LMIs
LMIs = [LMIs, W+W'>=theta2*eye(nx)];
for i=1:r
    SizeConstr = [theta1*eye(nx), Y{i}';
        Y{i},    eye(nu)];
    LMIs = [LMIs, SizeConstr>=0];
    for j=1:r
        if(i==j)
            continue;
        end
        for l=1:eta
            if(Modeltype == 2) %model with two indices
                LMIs = [LMIs, Upsilon{i,i,l}<=0];
                LMIs = [LMIs, 2/(r-1)*Upsilon{i,i,l} + ...
                    Upsilon{i,j,l} + Upsilon{j,i,l}<=0];
            else
                LMIs = [LMIs, Upsilon{i,l}<=0];
            end

        end
    end
end

% Configurando o Solver.
opts=sdpsettings;
% opts.solver='lmilab';
opts.solver='sedumi';
opts.verbose=0;
% Resolvendo as LMIs
sol = solvesdp(LMIs,[],opts);
che=min(checkset(LMIs));
if che > 0
    % Encontra o valor num�rico das matrizes
    W=double(W);
    K = cell(1,r);
    for i=1:r
        K{i} = double(Y{i})/W;
    end
else
    error('LMIs infact�veis')
    Q=[]; W=[]; K=[];
end
end
%%
function G=PdotConvex(p)

% Local models
ri=2^p;

% Total of columns possibilities
m2=2^ri;

% Number of columns of Total matrix
co=factorial(2^p)/(factorial(2^p/2)^2);


th=0:1:m2-1;

Vaux1=cell(1,m2);

for i=1:length(th)
    Vaux1{1,i}=de2bi([th(i)],ri);
end

Vaux2=zeros(ri,m2);
for j=1:m2
    for i=1:ri
        Vaux2(i,j)=(-1)^Vaux1{1,j}(i);
    end
end

% Total matrix
G=[];

for j=1:m2
    if sum(Vaux2(:,j))==0
        G=[G, Vaux2(:,j)];
    end
end
end

%%
function [K,P,R,L,A,G,Rset] = CNMAC2023(Modeltype,A,B)
if Modeltype ~= 1
    K = []
else
    Rset = 1:size(B,2);
    n = size(A{1},2);
    G=[1,2];
    % K calculated using Mozelli
    K{1} = [
    9.9990   -0.2241   -0.0000    0.0000    1.5038    0.0001    0.0000    0.0000
   -0.2241    9.7749   -0.0000    0.0000   -0.0001    1.5037   -0.0000    0.0000
    0.0000   -0.0000   -1.2490    0.0000    0.0000   -0.0000    0.9293    0.0000
   -0.0000   -0.0000   -0.0000   -1.4983   -0.0000   -0.0000   -0.0000    0.9294]
    K{2} = [
    9.9989    0.2241   -0.0000    0.0000    1.5032    0.0001    0.0000    0.0000
    0.2241    9.7748   -0.0000    0.0000   -0.0001    1.5031   -0.0000    0.0000
    0.0000    0.0000   -1.2602    0.0000    0.0000    0.0000    0.9281    0.0000
   -0.0000   -0.0000   -0.0000   -1.5095   -0.0000   -0.0000   -0.0000    0.9282]
   K{3} = [
    9.7747   -0.2241   -0.0000    0.0000    1.5024    0.0001    0.0000    0.0000
   -0.2241    9.9988    0.0000   -0.0000   -0.0001    1.5023   -0.0000    0.0000
    0.0000   -0.0000   -1.2757    0.0000    0.0000   -0.0000    0.9263    0.0000
   -0.0000   -0.0000   -0.0000   -1.5250   -0.0000   -0.0000   -0.0000    0.9265]
   K{4} = [
    9.7746    0.2241   -0.0000    0.0000    1.5009    0.0001    0.0000    0.0000
    0.2241    9.9987   -0.0000   -0.0000   -0.0001    1.5007   -0.0000    0.0000
    0.0000    0.0000   -1.2991    0.0000    0.0000    0.0000    0.9234    0.0000
   -0.0000   -0.0000   -0.0000   -1.5483   -0.0000   -0.0000   -0.0000    0.9236]
    for j = Rset
        A{j} = A{j}+B{j}*K{j};
    end
    h{1} = @(psi) (sin(2*psi)/4 - 1/2)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 1285/747);
    h{2} = @(psi) -(sin(2*psi)/4 + 1/2)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 1285/747);
    h{3} = @(psi) -(sin(2*psi)/4 - 1/2)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 538/747);
    h{4} = @(psi) (sin(2*psi)/4 + 1/2)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 538/747);

    dh{1} = @(psi) [ 0, 0, 0, 0, 0, 0, 0,(cos(2*psi)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 1285/747))/2 + 2*cos(psi)*sin(psi)*(sin(2*psi)/4 - 1/2)];
    dh{2} = @(psi) [ 0, 0, 0, 0, 0, 0, 0, - (cos(2*psi)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 1285/747))/2 - 2*cos(psi)*sin(psi)*(sin(2*psi)/4 + 1/2)];
    dh{3} = @(psi) [  0, 0, 0, 0, 0, 0, 0,- (cos(2*psi)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 538/747))/2 - 2*cos(psi)*sin(psi)*(sin(2*psi)/4 - 1/2)];
    dh{4} = @(psi) [ 0, 0, 0, 0, 0, 0, 0, (cos(2*psi)*((538*cos(psi)^2)/747 + (1285*sin(psi)^2)/747 - 538/747))/2 + 2*cos(psi)*sin(psi)*(sin(2*psi)/4 + 1/2)];
 
   

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

%     % Set estimation
%     V = @(x1,x2,x3,x4,x5,x6,x7,psi) sum(arrayfun(@(k) [x1,x2,x3,x4,x5,x6,x7,psi]*h{k}(psi)*P{k}*[x1,x2,x3,x4,x5,x6,x7,psi]',G));
%     hdot = @(x1,x2,x3,x4,x5,x6,x7,psi,k) sum(arrayfun(@(j) dh{k}(psi)*h{j}(psi)*A{j}*[x1,x2,x3,x4,x5,x6,x7,psi]',Rset));
%     Dset = @(x1,x2,x3,x4,x5,x6,x7,psi) sum(arrayfun(@(k) [x1,x2,x3,x4,x5,x5,x7,psi]*hdot(x1,x2,x3,x4,x5,x6,x7,psi,k)*P{k}*[x1,x2,x3,x4,x5,x6,x7,psi]',G));
% 
%     %calculate V and D
%     meshPoints=500;
%     tol=10/meshPoints;
%     x1 = linspace(-5,5,meshPoints);
%     x2 = linspace(-5,5,meshPoints);
%     x3 = linspace(-5,5,meshPoints);
%     x4 = linspace(-5,5,meshPoints);
%     x5 = linspace(-5,5,meshPoints);
%     x6 = linspace(-5,5,meshPoints);
%     x7 = linspace(-5,5,meshPoints);
%     psi = linspace(-pi,pi,meshPoints);
%     %Ve=zeros(length(x1),length(x2),length(x3),length(x4),length(x5),length(x6),length(x7),length(psi));
%     %De=zeros(length(x1),length(x2),length(x3),length(x4),length(x5),length(x6),length(x7),length(psi));
%     for X1=1:length(x1)
%         for X2=1:length(x2)
%             for X3=1:length(x3)
%                 for X4=1:length(x4)
%                     for X5=1:length(x5)
%                         for X6=1:length(x6)
%                             for X7=1:length(x7)
%                                 for X8=1:length(psi)
%                                     Ve(X1,X2,X3,X4,X5,X6,X7,X8) = V(x1(X1),x2(X2),x3(X3),x4(X4),x5(X5),x6(X6),x7(X7),psi(X8));
%                                     De(X1,X2,X3,X4,X5,X6,X7,X8) = Dset(x1(X1),x2(X2),x3(X3),x4(X4),x5(X5),x6(X6),x7(X7),psi(X8));
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
% 
%     for  comb = nchoosek(1:8,2)
%         
% 
%     end
% 
%     meshPoints=500;
%     tol=10/meshPoints;
%     x = linspace(-5,5,meshPoints);
%     y = linspace(-5,5,meshPoints);
% 
%     [X,Y]=meshgrid(x,y);
%     for i=1:length(x)
%         for j = 1:length(y)
%             Ve(i,j) = V(X(i,j),Y(i,j));
%             De(i,j) = Dset(X(i,j),Y(i,j));
%         end
%     end
% 
%     %calculate b
%     b=min([min(Ve(:,1)), min(Ve(:,end)), min(Ve(1,:)), min(Ve(end,:))])
%     b=fix(b*1e2)/1e2;
% 
%     figure(1);
%     [~,c]=contour(X,Y,Ve,linspace(0,b,5),'r','ShowText','on','DisplayName','V')
%     hold on
%     [~,d]=contour(X,Y,De,[0,fix(max(max(De))*1e2)/1e2],'b','ShowText','on','DisplayName','D')
%     legend;
%     %from the graph
% %     l = 0.18;
%     K = K;
end
end

