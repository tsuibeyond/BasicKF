clear
clc;

%%%%%%%%%%%%%%% 间接法 卡尔曼滤波 （正弦波跟踪）%%%%%%%%%%%%%%%%%%%%%%
t=0:0.01:(10-0.01);
kf.trueValue = 10*sin(0.5*pi*t);                      % 生成正弦波数据
kf.X_ins = zeros(size(kf.trueValue));
kf.X_err = zeros(size(kf.trueValue));
kf.measurement = 0.5*randn(size(t)) + kf.trueValue;     % 生成测量数据

nStates = 1;
mMeasures = 1;
kf.X_ins(1) = 0;   % 状态初始值为 0
kf.X_err(1) = 0;
kf.P = zeros(nStates,nStates,length(t));
kf.P(:,:,1) = 10;     % 状态估计协方差初始值
kf.Qt = 2.5;        % 模型噪声
kf.R = 0.5;         % 观测噪声
kf.K = zeros(nStates,mMeasures,length(t));  % 滤波增益

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
    
    % 一步预测 
%     kf.X_ins(i) = kf.PHI*kf.X_ins(i-1);
    kf.X_ins(i) = kf.X_ins(i-1) + kf.xdot*det_t;
    kf.X_err(i) = 0;
    kf.P(:,:,i) = kf.PHI*kf.P(:,:,i-1)*kf.PHI' + kf.Qk; % 一步预测协方差
    kf.K(:,:,i) = kf.P(:,:,i)*kf.H'/(kf.H*kf.P(:,:,i)*kf.H' + kf.R); % 滤波增益计算
    kf.Z = kf.measurement(i) - kf.H*kf.X_ins(i); 
    kf.X_err(i) = kf.X_err(i) + kf.K(:,:,i) * (kf.Z - kf.X_err(i));
    kf.P(:,:,i) = (eye(nStates) - kf.K(:,:,i)*kf.H) * kf.P(:,:,i);
    
    kf.X_ins(i) = kf.X_ins(i) + kf.X_err(i);
end

%%%%%%%%%%%滤波结果%%%%%%%%%%%%%%%%%%%%%%%%
smooth_res = medfilt1(kf.measurement,10); % 10点滑动平均滤波器（）
figure(1);
plot(t,kf.trueValue,'r',t,kf.X_ins,'g',t,kf.measurement,'b',t,smooth_res,'k');
legend('true value','kf estimation','measurement','10 point mean-filter');
xlabel('Sample time');
ylabel('Value');
title('indirect kalman filter');

rms_indirect = std(kf.X_ins-kf.trueValue)