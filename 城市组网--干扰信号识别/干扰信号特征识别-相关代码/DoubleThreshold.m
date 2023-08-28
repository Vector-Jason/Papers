
clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%@para
%
% jammerType
%
%       singleTone            ->        1     
%       multiTone             ->        2
%       linear sweep          ->        3
%       AM                    ->        4 
%       FM                    ->        5
%       NB AWGN               ->        6  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����ź�
jammerType = 6;            
jammerSignals = jammerSigFunc(jammerType);
jammerSignals = jammerSignals./abs(max(jammerSignals));
max(jammerSignals)

% jammerSignals = jammerSignals(2000:7000);

%% �ŵ�
JNR=-5;
cap_success = zeros(1,length(JNR));
for snr_index = 1:length(JNR)
for carlo_index = 1:1

rsignal=awgn(jammerSignals,JNR(snr_index),'measured');
L=length(rsignal);

j_f = abs(fft(jammerSignals));
j_f = j_f(1:L/2);
Rtemp=abs(fft(rsignal));
Rtemp=Rtemp./mean(Rtemp);
R=Rtemp(1:L/2);


%% Ԥ����
L_window=10;
R_smooth=window_smooth(R,L_window);

figure(22)

noise = rsignal - jammerSignals;
plot(1:length(noise),abs(noise));
hold on
plot(1:length(jammerSignals),abs(jammerSignals));
legend('Noise signal','Send signal');
ylabel('ʱ���ֵ');xlabel('sampling index');title(['JNR = ',num2str(JNR),' dB ']);

%% ˫�����������

% R_smooth = R_smooth.^0.5;
[R_sorted,order]=sort(R_smooth);%order�Ƕ�Ӧԭ���е�˳��
figure(2)
bar(1:length(R_sorted),R_sorted);
V_spectrum=[R_sorted;order];
Q=0.01;%��Ϊ��С��һ����������û�и��ŵ�
S_nointerf=sum(R_sorted(1:Q*length(R_sorted)));
length_S=length(R_sorted(1:Q*length(R_sorted)));

p_low=10e-5;%�������龯����
p_high=10e-8;%�������龯����
T_low=(2/sqrt(pi))*sqrt(-log10(p_low));
T_high=(2/sqrt(pi))*sqrt(-log10(p_high));%�������Ӽ���
% T_low = erfc(p_low);
% T_high = erfc(p_high);

V_remain=V_spectrum(:,Q*length(R_sorted)+1:end);%��Ҫ���Ƚϵ����߼���
S_nointerf_temp=0;
while 1
    E_nointerf=S_nointerf/length_S;
    THlow=E_nointerf*T_low;
    S_nointerf_temp=S_nointerf;
    for index=1:length(V_remain(1,:))%��ѭ������ѡ��������������ֵ������
        if(V_remain(1,index)<=THlow)
            S_nointerf=S_nointerf+V_remain(1,index);
            length_S=length_S+1;
            V_remain(1,index)=-1;%С������ֵ�����ߣ���ֵ-1�Ա�ǣ������޳�
        end
    end
    
    V_remain_next=zeros(2,length(V_remain(1,:)));
    index_temp=1;
    for index=1:length(V_remain(1,:))%��ѭ�������޳��ϸ�ѭ��ѡ��������
        if(V_remain(1,index)>0)
            V_remain_next(:,index_temp)=V_remain(:,index);
            index_temp=index_temp+1;
        end
    end
    %���´�������ֵ�����߼���
    if(index_temp>1)
        V_remain=V_remain_next(:,1:index_temp-1);
    else
        V_remain=[0;0];
    end
    %�������޸������ߣ�����ѭ��
    if ((S_nointerf_temp-S_nointerf))/S_nointerf_temp < 0
        break ;
    end
end


%�������ж�
    THhigh=E_nointerf*T_high;
    for index=1:length(V_remain(1,:))
        if(V_remain(1,index)<THhigh)
            V_remain(1,index)=-1;
        end
    end
    
    V_remain_next=zeros(2,length(V_remain(1,:)));index_temp=1;
    for index=1:length(V_remain(1,:))
        if(V_remain(1,index)>0)
            V_remain_next(:,index_temp)=V_remain(:,index);
            index_temp=index_temp+1;
        end
    end
    %����ʣ�����߼���
    if(index_temp>1)
        V_remain=V_remain_next(:,1:index_temp-1);
    else
        V_remain=[0;0];
    end


if V_remain ~= [0;0]
    cap_success(snr_index) = cap_success(snr_index) + 1;
else
    cap_success(snr_index) = cap_success(snr_index);
end
    

end
end

figure(3)
subplot(211)
semilogy(1:length(j_f(1:L/2)),abs(j_f(1:L/2)));axis([0 length(R_smooth)-1 10e-5 max(abs(j_f))*2]);
legend('Send signal');
ylabel('Ƶ���ֵ');xlabel('frequence index');title(['JNR = ',num2str(JNR),' dB']);

subplot(212)
semilogy(1:length(R_smooth),R_smooth);axis([0 length(R_smooth)-1 10e-5 max(R_smooth)*2]);
ylabel('Ƶ���ֵ');xlabel('frequence index');title(['JNR = ',num2str(JNR),' dB     Lwindow = ',num2str(L_window)]);
hold on
line([0 length(R_smooth)-1],[THlow THlow],'Color','red','LineStyle','--');
hold on
line([0 length(R_smooth)-1],[THhigh THhigh],'Color','cyan','LineStyle','--');
scatter(V_remain(2,:)-1,V_remain(1,:));
legend('Received signal','THlow','THhigh','Deteced signal','location','SouthEast');

