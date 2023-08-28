clear
clc
close all

% jammerType
%
%       singleTone            ->        1    单音 
%       multiTone             ->        2    多音
%       linear sweep          ->        3    线性扫频
%       AM                    ->        4    噪声调幅
%       FM                    ->        5    噪声调频
%       NB AWGN               ->        6    窄带高斯

jammerType  = 2;            
jammerSignals = jammerSigFunc(jammerType);

Y = awgn(jammerSignals,20,'measured');     %模拟awgn信道，加入噪声，达到不同的干噪比，Y为接收到的信号
Y = Y/max(Y);
% Y = abs(Y);

F = fft(Y); % 接收到的信号进行傅里叶变换
F = abs(F);
F = F/max(F); % 归一化jammerSignals

%% -------------时域 特征提取-------------------%%
pd = skewness(Y) % 时域矩偏度：信号的分布偏离中心的程度；大于0说明信号在右侧的部分多于左侧，反之小于0
fd = kurtosis(Y); % 时域矩峰度：信号分布的尖锐程度；大于3说明波形出现多个尖峰，整体分布较为陡峭；小于3说明信号分布较为平缓，呈现扁平分布
mea = mean(Y(:)); % 均值
fc = var(Y); % 方差
R = fc/(mea.^2); % 时域包络起伏度：信号振幅随时间变化的剧烈程度

%% -------------时频域特征提取------------------%%
Rf = zeros(1,20000); max1 = zeros(5);
for a=0:0.5:2
    b=a*2+1;
    Rf = myfrft(Y,a);
    max1(b)=abs(max(Rf));
end
M = max(max1,[],'all'); % 分数阶傅里叶域最大值

%% -------------波形域特征提取------------------%%
V = Y;
N = length(V);
V(N+1) = 0; d2 = 0; d1 = 0;
for i = 1:N
    d1 = d1+abs(V(i)-V(i+1));
end
for i = 1:(N/2)
    max2 = max(max(V(2*i-1),V(2*i)),V(2*i+1));
    min2 = min(min(V(2*i-1),V(2*i)),V(2*i+1));
    d2 = d2+(max2-min2);
    max2 = 0; min2 = 0;
end
Df = 1 + (log(d1/d2))/(log(2)); % 盒维数：描述分形信号的几何尺度信息

%% -------------频域 特征提取------------------%%
pd2 = skewness(F); % 频域矩偏度：对正态分布的偏离程度，信号频域的对称性
fd2 = kurtosis(F); % 频域矩峰度：信号的陡峭程度
crestfactor=max(F)/mean(F); % 频谱峰度：信号频谱的陡峭程度
[max, m]=max(F); a = 0; b = 0; d = 0;

% -------------计算单频能量聚集度--------------------
for i=1:N
    a = a + F(i)^2;
end
C = (F(m)^2+F(m+1)^2) / a; % 单频能量聚集度

% -------------计算平均频谱平坦系数------------------
F1 = zeros(1,N); % 分配空间
for k = 1:N
    %计算和式
    if k <= 600
    temp_sum = sum(F(1:600+k));
    elseif k >= 19401
    temp_sum = sum(F(k-600:N));
    else
    temp_sum = sum(F(k-600:k+600));
    end
    %计算F1(k)
    F1(k) = F(k) - temp_sum/1201; 
end
%平均频谱平坦系数：表现信号频谱幅度的变化程度
Fc = sqrt(sum((F1-mean(F1)).^2)/N);