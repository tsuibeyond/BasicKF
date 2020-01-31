% clear
% clc;

%%%%%%%%%%%%%%%ֱ�ӷ� �������˲� �����Ҳ����٣�%%%%%%%%%%%%%%%%%%%%%%
t=0:0.01:(10-0.01);
kf.trueValue = 10*sin(0.5*pi*t);             % �������Ҳ�����
kf.X = zeros(size(kf.trueValue));
kf.Z = kf.measurement; %0.5*randn(size(t)) + kf.trueValue;    % ���ɲ�������

nStates = 1;
mMeasures = 1;
kf.X(1) = 0;   % ״̬��ʼֵΪ 0
kf.P = zeros(nStates,nStates,length(t));
kf.P(:,:,1) = 10;     % ״̬����Э�����ʼֵ
kf.Qt = 2.5;    % ģ������
kf.R = 0.5;    % �۲�����
kf.K = zeros(nStates,mMeasures,length(t));  % �˲�����

det_t = t(2)-t(1);
for i = 2 : length(t)
    kf.xdot = 0;
    kf.F = 0;
    kf.H = 1;
    AA = [-kf.F kf.Qt; 
        zeros(nStates,nStates) kf.F']*det_t;
    BB = eye(nStates*2) + AA;
    kf.PHI = BB(nStates+1:2*nStates,nStates+1:2*nStates)';
    kf.Qk = kf.PHI*BB(1:nStates,nStates+1:2*nStates);
    
%     kf.X(i) = kf.PHI*kf.X(i - 1);   % һ��Ԥ�� 
    kf.X(i) = kf.X(i-1) + kf.xdot*det_t;
    kf.P(:,:,i) = kf.PHI*kf.P(:,:,i-1)*kf.PHI' + kf.Qk; % һ��Ԥ��Э����
    kf.K(:,:,i) = kf.P(:,:,i)*kf.H'/(kf.H*kf.P(:,:,i)*kf.H' + kf.R); % �˲��������
    kf.X(i) = kf.X(i) + kf.K(:,:,i) * (kf.Z(i) - kf.X(i));
    kf.P(:,:,i) = (eye(nStates) - kf.K(:,:,i)*kf.H) * kf.P(:,:,i);
end

%%%%%%%%%%%�˲����%%%%%%%%%%%%%%%%%%%%%%%%
smooth_res = medfilt1(kf.Z,10); % 10�㻬��ƽ���˲�������
figure(1);
plot(t,kf.trueValue,'r',t,kf.X,'g',t,kf.Z,'b',t,smooth_res,'k');
legend('true value','kf estimation','measurement','10 point mean-filter');
xlabel('Sample time');
ylabel('Value');
title('direct kalman filter');

rms_direct = std(kf.X-kf.trueValue)

