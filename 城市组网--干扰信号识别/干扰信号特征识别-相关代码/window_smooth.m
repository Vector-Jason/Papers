function R_smooth = window_smooth(R,L)
%����������������������������������������
%WINDOW_SMOOTH ������ֵ�˲�
%RΪ��Ҫ���������У�LΪ�������ڵĳ���
%����������������������������������������
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
R_smooth=[Rtemp R_smooth(L:length(R))];%���е�ǰ�������ǻ�������δ��ȫ�������У�Ĭ��ΪǰL������ƽ��

% figure('Name','window_smooth')
% plot(0:length(R_smooth)-1,10*log10(R_smooth));axis([0 length(R_smooth) -40 max(10*log10(R_smooth))*2]);
% title('����ƽ���˲���ķ�����');ylabel('���ʣ�dBW��');xlabel('Ƶ�ʣ�Hz��');

end

