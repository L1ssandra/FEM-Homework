%hw2auto.m  �Զ�����������

%% �������
L = 4; % �ܹ��ļ������
Table = zeros(L,4); % ���

N0 = 10; % ��ʼ���ʷֶ���

%% ���
for k = 1:L
    
    N = N0*2^(k - 1); % ��������
    
    hw2 % ������
    
    Table(k,1) = L2; % ��L2�������
    Table(k,3) = H1; % ��H1�������
    
    if k > 1 % ȷ���о�����
        Table(k,2) = log(Table(k - 1,1)./Table(k,1))./log(2); % ����L2 Error�Ľ�
        Table(k,4) = log(Table(k - 1,3)./Table(k,3))./log(2); % ����H1 Error�Ľ�
    end
    
end