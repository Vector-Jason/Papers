

function jammerSignals = jammerSigFunc(jammerType)

    switch jammerType
        
        case 1
                    Fs=20000;  %采样频率
                    N=20000;    %采样点
                    n=0:N-1;
                    t=n/Fs;  %时间序列
                    fc=1000;  %载波信号频率 
                    f=n*Fs/N;  %频率 
                    Uc=1*exp(2*1i*fc*pi*t);     %载波信号 
                    C1=fft(Uc);             %对载波信号进行傅里叶变换 
                    cxf=abs(C1);           %进行傅里叶变换  
                    cxf=cxf/max(cxf);%归一化
                    subplot(2,1,1);plot(t,real(Uc));title('载波信号波形');xlabel('时间(s)');ylabel('幅度(V)');title('单音干扰信号波形');axis([0 0.009 -1 1]);
                    subplot(2,1,2); plot(f(1:N/2),cxf(1:N/2));title('载波信号频谱'); axis([0 2000 0 1]);xlabel('频率(Hz)');ylabel('功率');title('单音干扰信号归一化功率谱');
                    set(gca,'YTick',0:1:1);
                    jammerSignals = Uc;

        case 2
                    Fs=20000;  
                    N=20000;    
                    n=0:N-1;
                    t=n/Fs;  
                    A0=1;  
                    fc=1000;  
                    f=n*Fs/N; 
                    w0=2*fc*pi*1i; 
                    step=2*pi*50*1i;
                    Uc=A0*exp(w0*t)+A0*exp((w0+step)*t)++A0*exp((w0+2*step)*t)++A0*exp((w0+3*step)*t)+A0*exp((w0-step)*t)++A0*exp((w0-2*step)*t)++A0*exp((w0-3*step)*t);%多音信号 
                    C1=fft(Uc);      
                    cxf=abs(C1);     
                    cxf=cxf/max(cxf);
                    subplot(2,1,1);
                    plot(t,real(Uc));
                    xlabel('时间(s)');
                    ylabel('幅度(V)');
                    title('多音干扰波形');
                    axis([0 0.1 -8 8]);
                    subplot(2,1,2);
                    plot(f(1:N/2),cxf(1:N/2));
                    title('多音干扰频谱');
                    axis([0 2000 0 1]);
                    xlabel('频率(Hz)');
                    ylabel('功率');
                    title('多音干扰归一化功率谱');
                    set(gca,'YTick',0:1:1);
                    jammerSignals = Uc;

        case 3

                    t=0:0.0001:1-0.0001;
                    f0=300;%扫频起始频率
                    fe=2000;%扫频截止频率
                    x=chirp(mod(t,1),f0,1,fe);
                    subplot(3,1,1);plot(t,x);title('线性扫频信号波形');xlabel('时间(s)');ylabel('幅度(V)');
                     
                    ft=f0+(fe-f0)*mod(t,1);
                    subplot(3,1,2);plot(t,ft);title('线性扫频信号频率-时间图');xlabel('时间(s)');ylabel('频率(Hz)');
                     
                    t=0:0.0001:1-0.0001;
                    x=chirp(t,f0,1,fe);
                    C1=fft(x);     
                    cxf=abs(C1);    
                    cxf=cxf/max(cxf);
                    subplot(3,1,3);plot(cxf); axis([0 length(t)-1 0 1]);title('线性扫频信号归一化频谱');xlabel('频率(Hz)');ylabel('功率');
                    jammerSignals = x;

        case 4

                    fj=20e6;
                    fs=4*fj; 
                    Tr=200e-6;
                    
                    t1=0:1/fs:1.25*Tr-1/fs; 
                    N=length(t1);
                    u=wgn(1,N,0);
                    
                    df1=fs/N;n=0:N/2;f=n*df1;
                    
                    wp=10e6;
                    ws=14e6;
                    rp=1; 
                    rs=60;
                    [n1,wn1]=buttord(wp/(fs/2),ws/(fs/2),rp,rs);
                    [b,a]=butter(n1,wn1);
                    u1=filter(b,a,u);
                    
                    p=0.1503*mean((u1.^2));
                    figure;
                    subplot(2,2,1),
                    plot(t1,u1),
                    title('噪声信号波形'); axis([0,0.02e-4,-2,2]);xlabel('时间(s)');ylabel('幅度(V)');
                    subplot(2,2,2)
                    j2=fft(u1);
                    plot(f,10*log10(abs(j2(n+1)*2/N)));
                    xlabel('频率(Hz)');ylabel('功率(dBW)');axis([0,4e7,-70,0]);title( '噪声信号功率谱');
                    
                    u0=1;
                    y=(u1+u0).*cos(2*pi*fj*t1+2);
                    u2=u1+u0;
                    u3=-u0-u1;
                    subplot(2,2,3), plot(t1,y,t1,u2,t1,u3),title( '噪声调幅信号时域波形'); axis([0,0.02e-4,-2,2]);xlabel('时间(s)');ylabel('幅度(V)');
                    subplot(2,2,4), J=fft(u1);plot(f,10*log10(abs(J(n+1))));xlabel('频率(Hz)');ylabel('功率(dBW)');axis([0,4e7,-20,50]);title( '噪声调幅信号功率谱');
                    jammerSignals = u1;
            
            
        case 5
                    uj=1;mf=2;wpp=10;
                    
                    fj=20e6;
                    fs=8*fj;
                    Tr=520e-6;
                    t1=0:1/fs:3*Tr-1/fs;
                    N=length(t1);
                    u=wgn(1,N,0);
                    wp=10e6;
                    ws=16e6;
                    rp=1;
                    rs=60;
                    [n1,wn1]=buttord(wp/(fs/2),ws/(fs/2),rp,rs);
                    [b,a]=butter(n1,wn1);
                    u1=filter(b,a,u);
                    p=0.8503*mean((u1.^2)) ;

                    
                    fj=20e6;
                    fs=8*fj;
                    Tr=520e-6;
                    bj=5e6;
                    t1=0:1/fs:3*Tr-1/fs;N=length(t1);
                    u=wgn(1,N,wpp);
                    df1=fs/N;
                    n=0:N/2;
                    f=n*df1;
                    wp=10e6;
                    ws=14e6;
                    rp=1;
                    rs=60;
                    [Nn,wn]=buttord(wp/(30e6/2),ws/(30e6/2),rp,rs);
                    [b,a]=butter(Nn,wn);


                    figure;subplot(2,2,1),plot(t1,u1),title('噪声信号波形');axis([0,2e-6,-2,2]);xlabel('时间(s)');ylabel('幅度(V)');
                    subplot(2,2,2),j2=fft(u1); plot(f,10*log10(abs(j2(n+1)*2/N)));xlabel('频率(Hz)');ylabel('功率（dBW）');axis([0,4e7,-20,50]);title( '噪声信号功率谱');axis([0,4e7,-80,0]);
                    i=1:N-1;
                    ss=cumsum([0 u1(i)]);
                    ss=ss*Tr/N;
                    y=uj*cos(2*pi*fj*t1+2*pi*mf*bj*ss*10);
                    subplot(2,2,3), plot(t1,y),title( '噪声调频信号波形'),axis([0,2e-6,-1.5,1.5]);xlabel('时间(s)');ylabel('幅度(V)');
                    y=uj*cos(2*pi*fj*t1+2*pi*mf*bj*ss);  
                    subplot(2,2,4),J=fft(y);plot(f,10*log10(abs(J(n+1))));axis([0,4e7,-20,60]);xlabel('频率(Hz)');ylabel('功率（dBW）');axis([0,4e7,-20,50]);title( '噪声调频信号功率谱')
                    jammerSignals = y(1:20000);
%                     for index = 1:20000
%                          jammerSignals(index) = y(index*12) ;
% 
%                     end


        case 6

                    fs=1000;%采样频率hz
                    T_N=20;%总时间s
                    t=1/fs:1/fs:T_N;%时间向量
                    L=T_N*fs;%样本数量
                    power=3;%噪声功率,单位为dbw
                    z=wgn(L,1,power);
                    figure(1)
                    subplot(2,2,1)
                    plot(t,z)
                    xlabel("时间/s")
                    ylabel("幅度/v")
                    title("高斯白噪声（时域）")
                    
                    fft_z=fft(z);%快速傅里叶变换之后的噪声
                    P = abs(fft_z/L);%取幅频特性，除以L
                    P = P(1:L/2+1);%截取前半段
                    P(2:end-1)=2*P(2:end-1);%单侧频谱非直流分量记得乘以2
                    f = fs*(0:(L/2))/L;%奈奎斯特采样定理
                    subplot(2,2,2)
                    plot(f,P)
                    xlabel("频率/Hz")
                    ylabel("幅度/v")
                    title("高斯白噪声（频域）")
                    
                    [b,a]=butter(8,[300/(fs/2),400/(fs/2) ]);%获得8阶巴特沃斯滤波器系数，100-200Hz

%                     freqs(b,a)%画滤波器特性曲线
                    lvbo_z=filter(b,a,z);%滤波
                    

                    subplot(2,2,3)
                    plot((lvbo_z))
                    xlabel("时间/Hz")
                    ylabel("幅度/v")
                    title("窄带高斯噪声（时域）")
                    
                    fft_lvbo_z=fft(lvbo_z);%傅里叶变换
                    P = abs(fft_lvbo_z/L);%取幅频特性，除以L
                    P = P(1:L/2+1);%截取前半段
                    P(2:end-1)=2*P(2:end-1);%单侧频谱非直流分量记得乘以2
                    subplot(2,2,4)
                    plot(f,P)
                    xlabel("频率/Hz")
                    ylabel("幅度/v")
                    title("窄带高斯噪声（频域）")
                    jammerSignals = lvbo_z.';

    end
            
        


end