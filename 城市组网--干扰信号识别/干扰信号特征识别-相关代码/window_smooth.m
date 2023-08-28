function R_smooth = window_smooth(R,L)
%————————————————————
%WINDOW_SMOOTH 滑动均值滤波
%R为需要操作的序列，L为滑动窗口的长度
%————————————————————
c=zeros(1,L);
for index=1:L
    c(index)=1/L;
end
R_smooth=conv(c,R);
% R_smooth=R_smooth(1:length(R));
Rtemp=zeros(1,L-1);
for index=1:length(Rtemp)
    Rtemp(index)=R_smooth(L);
end
R_smooth=[Rtemp R_smooth(L:length(R))];%序列的前几个数是滑动窗口未完全进入序列，默认为前L个数的平均

% figure('Name','window_smooth')
% plot(0:length(R_smooth)-1,10*log10(R_smooth));axis([0 length(R_smooth) -40 max(10*log10(R_smooth))*2]);
% title('滑动平均滤波后的幅度谱');ylabel('功率（dBW）');xlabel('频率（Hz）');

end

