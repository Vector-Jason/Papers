function Faf = myfrft(f, a)
%分数阶傅里叶变换函数
%输入参数f为原始信号，a为阶数
%输出结果为原始信号的a阶傅里叶变换
N = length(f);%总采样点数
shft = rem((0:N-1)+fix(N/2),N)+1;%此项等同于fftshift(1:N)，起到翻转坐标轴的作用
sN = sqrt(N);%对离散傅里叶变换的定义，有这个乘积项
a = mod(a,4);%按周期性将a定义在[0,4]

%特殊情况直接处理
if (a==0), Faf = f; return; end%自身
if (a==2), Faf = flipud(f); return; end%f(-x)
if (a==1), Faf(shft,1) = fft(f(shft))/sN; return; end%f的傅里叶变换
if (a==3), Faf(shft,1) = ifft(f(shft))*sN; return; end%f的逆傅里叶变换

%利用叠加性将阶数变换到0.5 < a < 1.5
if (a>2.0), a = a-2; f = flipud(f); end%a=2是反转
if (a>1.5), a = a-1; f(shft,1) = fft(f(shft))/sN; end%a=1是傅里叶变换
if (a<0.5), a = a+1; f(shft,1) = ifft(f(shft))*sN; end%a=-1是逆傅里叶变换

%开始正式的变换
alpha = a*pi/2;
tana2 = tan(alpha/2);
sina = sin(alpha);

f = [zeros(N-1,1) ; interp(f) ; zeros(N-1,1)];%使用香农插值，拓展为4N
% 以下操作对应原论文中公式（29）
% 线性调频预调制
chrp = exp(-1i*pi/N*tana2/4*(-2*N+2:2*N-2)'.^2);
f = chrp.*f;
% 线性调频卷积
c = pi/N/sina/4;
Faf = fconv(exp(1i*c*(-(4*N-4):4*N-4)'.^2),f);
Faf = Faf(4*N-3:8*N-7)*sqrt(c/pi);
% 线性调频后调制
Faf = chrp.*Faf;
% 乘以最前面的A_Phi项
Faf = exp(-1i*(1-a)*pi/4)*Faf(N:2:end-N+1);

end

function xint=interp(x)%香农插值
% sinc interpolation
N = length(x);
y = zeros(2*N-1,1);
y(1:2:2*N-1) = x;
xint = fconv(y(1:2*N-1), sinc([-(2*N-3):(2*N-3)]'/2));%计算卷积
xint = xint(2*N-2:end-2*N+3);
end

function z = fconv(x,y)%利用fft快速计算卷积
N = length([x(:);y(:)])-1;%计算最大点数
P = 2^nextpow2(N);%补零
z = ifft( fft(x,P) .* fft(y,P));%频域相乘，时域卷积
z = z(1:N);%去零
end
