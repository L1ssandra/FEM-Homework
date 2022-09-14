% hw2.m  ������ű�

%% �������
xa = 0; % ������˵�
xb = 1; % �����Ҷ˵�
%N = 20; % �ʷֶ���������hw1autoʱע�͵�
h = (xb - xa)/N; % ���䳤��
X = xa:h:xb; % ����
XX = xa:h/2:xb; % ���е������
draw = 1; % ���Ϊ1�ͻ�ͼ

u = @(x) (x - 1).*sin(x); % ���
du = @(x) sin(x) + (x - 1).*cos(x); % ���ĵ���
f = @(x) -2.*cos(x) + (x - 1).*sin(x); % �Ҷ���

% [0,1]�ϵı�׼������
phi0 = @(x) (2.*x - 1).*(x - 1);
phi1 = @(x) 4.*x.*(1 - x);
phi2 = @(x) (2.*x - 1).*x;

% ����
dphi0 = @(x) 4.*x - 3;
dphi1 = @(x) 4 - 8.*x;
dphi2 = @(x) 4.*x - 1;

%% ��װ�Ҷ����� F = [ (f,phi_1/2),...,(f,phi_N - 1/2) ]'
F = zeros(2*N + 1,1);

for i = 1:N
    fp0 = @(x) f(h.*x + X(i)).*phi0(x).*h;
    fp1 = @(x) f(h.*x + X(i)).*phi1(x).*h;
    fp2 = @(x) f(h.*x + X(i)).*phi2(x).*h;
    b0 = quadgk(fp0,0,1);
    b1 = quadgk(fp1,0,1);
    b2 = quadgk(fp2,0,1);
    F(2*i - 1:2*i + 1) = F(2*i - 1:2*i + 1) + [b0;b1;b2]; % ��װ
end

F = F(2:end - 1); % ��ͷȥβ

%% ��װ�նȾ��� A = [ a(phi_i,phi_j) ]
% Ϊ�˼��ٷ�������⣬����ʹ��ϡ�����
A = sparse(2*N + 1,2*N + 1);

K = (1/(3*h))*sparse([7,-8,1;-8,16,-8;1,-8,7]); % �ֲ��նȾ���

for i = 1:N
    A(2*i - 1:2*i + 1,2*i - 1:2*i + 1) = A(2*i - 1:2*i + 1,2*i - 1:2*i + 1) + K; % ��װ
end

A = A(2:end - 1,2:end - 1); % ��ͷȥβ

%% �����Է�����
U = A\F;
U = [0;U;0];

%% ��ͼ
if draw == 1
    plot(XX,U,'r*',XX,u(XX),'b-')
end

%% �������
L2 = 0;
H1 = 0;
for i = 1:N
    % ����ת��Ϊ[0,1]�ϵĻ���
    e2 = @(x) h.*(U(2*i - 1).*phi0(x) + U(2*i).*phi1(x) + U(2*i + 1).*phi2(x) - u(h.*x + X(i))).^2;
    He2 = @(x) h.*(U(2*i - 1).*dphi0(x)/h + U(2*i).*dphi1(x)/h + U(2*i + 1).*dphi2(x)/h - du(h.*x + X(i))).^2 + e2(x);
    L2 = L2 + quadgk(e2,0,1);
    H1 = H1 + quadgk(He2,0,1);
end
L2 = sqrt(L2);
H1 = sqrt(H1);