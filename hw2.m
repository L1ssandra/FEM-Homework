% hw2.m  主计算脚本

%% 输入参数
xa = 0; % 区间左端点
xb = 1; % 区间右端点
%N = 20; % 剖分段数，调用hw1auto时注释掉
h = (xb - xa)/N; % 区间长度
X = xa:h:xb; % 网格
XX = xa:h/2:xb; % 加中点的网格
draw = 1; % 如果为1就画图

u = @(x) (x - 1).*sin(x); % 真解
du = @(x) sin(x) + (x - 1).*cos(x); % 真解的导数
f = @(x) -2.*cos(x) + (x - 1).*sin(x); % 右端项

% [0,1]上的标准基函数
phi0 = @(x) (2.*x - 1).*(x - 1);
phi1 = @(x) 4.*x.*(1 - x);
phi2 = @(x) (2.*x - 1).*x;

% 导数
dphi0 = @(x) 4.*x - 3;
dphi1 = @(x) 4 - 8.*x;
dphi2 = @(x) 4.*x - 1;

%% 组装右端向量 F = [ (f,phi_1/2),...,(f,phi_N - 1/2) ]'
F = zeros(2*N + 1,1);

for i = 1:N
    fp0 = @(x) f(h.*x + X(i)).*phi0(x).*h;
    fp1 = @(x) f(h.*x + X(i)).*phi1(x).*h;
    fp2 = @(x) f(h.*x + X(i)).*phi2(x).*h;
    b0 = quadgk(fp0,0,1);
    b1 = quadgk(fp1,0,1);
    b2 = quadgk(fp2,0,1);
    F(2*i - 1:2*i + 1) = F(2*i - 1:2*i + 1) + [b0;b1;b2]; % 安装
end

F = F(2:end - 1); % 掐头去尾

%% 组装刚度矩阵 A = [ a(phi_i,phi_j) ]
% 为了加速方程组求解，这里使用稀疏矩阵
A = sparse(2*N + 1,2*N + 1);

K = (1/(3*h))*sparse([7,-8,1;-8,16,-8;1,-8,7]); % 局部刚度矩阵

for i = 1:N
    A(2*i - 1:2*i + 1,2*i - 1:2*i + 1) = A(2*i - 1:2*i + 1,2*i - 1:2*i + 1) + K; % 安装
end

A = A(2:end - 1,2:end - 1); % 掐头去尾

%% 解线性方程组
U = A\F;
U = [0;U;0];

%% 画图
if draw == 1
    plot(XX,U,'r*',XX,u(XX),'b-')
end

%% 计算误差
L2 = 0;
H1 = 0;
for i = 1:N
    % 这里转化为[0,1]上的积分
    e2 = @(x) h.*(U(2*i - 1).*phi0(x) + U(2*i).*phi1(x) + U(2*i + 1).*phi2(x) - u(h.*x + X(i))).^2;
    He2 = @(x) h.*(U(2*i - 1).*dphi0(x)/h + U(2*i).*dphi1(x)/h + U(2*i + 1).*dphi2(x)/h - du(h.*x + X(i))).^2 + e2(x);
    L2 = L2 + quadgk(e2,0,1);
    H1 = H1 + quadgk(He2,0,1);
end
L2 = sqrt(L2);
H1 = sqrt(H1);