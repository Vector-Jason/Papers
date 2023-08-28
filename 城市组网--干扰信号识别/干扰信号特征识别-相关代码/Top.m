

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
%       partial-band          ->        4 
%       FM                    ->        5
%       NB AWGN               ->        6  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jammerTypeTemp = '单音干扰';'多音干扰';'线性扫频干扰';'窄带干扰';'噪声调频干扰';'噪声调幅干扰';

%% jammer signal generate

jammerType  = 6;            

jammerSignals = jammerSigFunc(jammerType);

figure(2)
% % tf-2D
% pspectrum(jammerSignals,2e3,'spectrogram','TimeResolution',0.1,'OverlapPercent',99,'Leakage',0.85)

% tf-3D
[p,f,t] =  pspectrum(jammerSignals,2e3,'spectrogram');
waterfall(f,t,p')
view([139,41])
xlabel('Frequency (Hz)')
ylabel('Time (seconds)')
% title(num2str(jammerTypeTemp(jammerType)))
wtf = gca;
wtf.XDir = 'reverse';



%% generate urben jamming channel




%% generate node received jamming signal



%% jamming feature extraction





















