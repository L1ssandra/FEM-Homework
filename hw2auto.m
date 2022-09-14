%hw2auto.m  自动计算收敛阶

%% 建立表格
L = 4; % 总共的计算次数
Table = zeros(L,4); % 表格

N0 = 10; % 初始的剖分段数

%% 填表
for k = 1:L
    
    N = N0*2^(k - 1); % 加密网格
    
    hw2 % 计算结果
    
    Table(k,1) = L2; % 把L2误差填入
    Table(k,3) = H1; % 把H1误差填入
    
    if k > 1 % 确保有旧数据
        Table(k,2) = log(Table(k - 1,1)./Table(k,1))./log(2); % 计算L2 Error的阶
        Table(k,4) = log(Table(k - 1,3)./Table(k,3))./log(2); % 计算H1 Error的阶
    end
    
end