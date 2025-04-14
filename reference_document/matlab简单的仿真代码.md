```matlab
%% initialization
clear; % 清除工作空间的所有变量
close all; % 关闭所有的Figure画布
clc; % 清除命令窗口的内容，对工作空间中的全部变量无任何影响

%% 创建一个输入信号(常数、比例函数、正弦函数、方波、锯齿波)
Fs = 1000; % 采样频率
t = 0 : 1/Fs : 1; % 时间步长,共1s长度
f = 50; % 信号频率为50Hz，w = 2*pi*f 机器人小陀螺转动频率约为(50Hz～100HZ)
k = 1; % 比例系数

constant_signal = ones(size(t)); % 创建一个常数信号，所有值都为1
ramp_signal = k*t;  % 创建一个线性递增的比例函数信号（对于simulink可用 'gain + clock'）
% 下列信号具有周期性质，故需要至少三个元素：频率f（周期T）、最大值Umax、最小值Umin
sine_signal = sin(2*pi*f*t) % 创建一个正弦信号(eg:sin(Wd*t) = sin(2*pi*f))
square_signal = square(2*pi*f*t);  % 创建一个方波信号
sawtooth_signal = sawtooth(2*pi*f*t); % 创建一个锯齿波信号

figure; % 创建画布(简而言之就是创建一个可以画图的背景版)
subplot(5,1,1); % 创建子图序列（无需多言）
plot(t, constant_signal) % 真正画图的函数
title('常数信号(t)') % 标题

subplot(5,1,2)
plot(t, ramp_signal)
title('斜坡信号(t)')

subplot(5,1,3);
plot(t, sine_signal);
title('正弦信号(t)');

subplot(5,1,4)
plot(t, square_signal);
title('方波信号(t)');

subplot(5,1,5)
plot(t, sawtooth_signal);
title('锯齿波信号(t)');

% ！傅里叶分析输入信号 ！
% 信号数组（上述已创建的信号集合）
signals = {constant_signal, ramp_signal, sine_signal, square_signal, sawtooth_signal}; % 元胞数组，取值一般使用{}花括号，()小括号是取地址赋值用的，可以把元胞数组类比成python中的列表
signal_names = {'Constant Signal', 'Ramp Signal', 'Sine Signal', 'Square Signal', 'Sawtooth Signal'};

% FFT分析和频谱图绘制(tips:套模板) 或者使用simulink中的PowerGui分析
for i = 1:length(signals) % 共计length个信号
    signal = signals{i}; % 取出信号的代表字符(i:待分析信号)
    figure;
    subplot(2,1,1);
    plot(t, signal); % 绘制信号
    title(signal_names{i});
    xlabel('Time (s)');
    ylabel('Amplitude');
    
    % FFT
    N = length(signal); % 信号长度
    Y = fft(signal); % FFT
    P2 = abs(Y/N); % 双边频谱
    P1 = P2(1:N/2+1); % 单边频谱
    P1(2:end-1) = 2*P1(2:end-1);
    
    % 傅里叶分析频率轴
    f_ = Fs*(0:(N/2))/N;
    
    % 绘制频谱图
    subplot(2,1,2);
    plot(f_, P1);
    title([signal_names{i} ' Fourier Spectrum']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
end

% tip:使用'gensig（）函数创建信号并使用[u,t]向量接收'——三要素：‘信号名字’，周期，x轴绘图长度，采样周期Fs(tups:离散化)
% 可以将MATLAB的M语言类比C++，与C++一样，M语言具有重载函数
T = 5;
T_Range = 30; % 显示（T_Range/T）个周期
Fs_gen = 0.01; % 离散采样周期(可省略)

gen_signals = {'sine', 'square', 'pulse'}
gen_signal_names = {'Sin Signal', 'Square Signal', 'Pulse Signal'};

for i = 1:length(gen_signals)
    sign =  gen_signals{i}
    [u,t] = gensig(sign, T, T_Range, Fs_gen);% 'sin'、'square'、'pulse',周期T,持续时间TLength，采样时间Fs
    figure;
    %subplot(length(gen_signals),1,1);
    plot(t, u); % 绘制信号
    title(gen_signal_names{i});
    xlabel('Time (s)');
    ylabel('Amplitude');
end
	
%% 创建系统传递函数并分析响应
% tf/zpk 连续与离散化
sys = tf([1 2],[1 3 5],0.01);
sys_zpk = zpk([1 2],[1 3 5],1,0.01)

% 延迟
sys_delay = tf([1 2],[1 3 5],'InputDelay',0.1); % tao = 0.1
% 其中:
sys_delay =
 
                    s + 2
  exp(-0.1*s) * -------------
                s^2 + 3 s + 5

% s 域创建
s = tf('s'); %创建关于传递函数的's'变量
Go = s+3/(s+1) ;%创建名为‘G’系统的传递函数
Gc = feedback(0.7*Go,1);
step(Go,Gc);
step(Gc);
impulse(Go);

% NUM/DEN 创建
num = [1 2];
den = [1 3 5];
sys = tf(num,den);

step(sys)
impulse(sys)

% 利用已知信号分析响应(以Go为例子)
input_signal = 5; % 选择输入信号
u = signals{input_signal}; % u 可以为任意函数与时间集合作用的输入集合
%[y,t] = lsim(Go,u,t);
lsim(Go,sys,u,t) % linear simuliation 传入传递函数、输入u 、时间t 返回[y,t]
xlabel('Time (s)');
ylabel('Amplitude');

%% 查看零极点分布
s = tf('s');
G = (s+1)/(s^2+s+1);
zpk_G = zpk(G)
pzmap(zpk_G);
[p,z] = pzmap(zpk_G);%获取极点零点
[z,p,k] = zpkdata(zpk_G,'v');

%% 查看零极点分布（以超前补偿器和滞后补偿器为例子）
Alpha = 1.5;
Beta = 1.5;
T = 1

s = tf('s');

Gc_Lead = (T*Alpha*s + 1)/(T*s + 1); % 超前补偿器
zpk_Gc_Lead = zpk(Gc_Lead)

pzmap(zpk_Gc_Lead);
[p,z] = pzmap(zpk_Gc_Lead);%获取极点零点
[z,p,k] = zpkdata(zpk_Gc_Lead,'v');% 获取极点零点增益

rlocus(Gc_Lead)
bode(Gc_Lead)

Gc_Lag = ((T*s + 1)/(Beta*T*s + 1))
bode(Gc_Lag)

%% PID控制器
s = tf('s');
Go = 1/s^2;
H = 1;
kp = 1;
ki = 2;
kd = 3;
Tf = 100
Gc_PID = pid(kp,ki,kd,Tf)
Gc = feedback(Go*Gc_PID,H)
bode(Gc_PID)


%% 分析根轨迹（开环->闭环）
rlocus(Go)
rlocus(sys)

%% 伯德图分析
bode(Go,Gc)
bode(sys)

%% 奈奎斯特分析
nyquist(Go)
nyquist(sys)

%% 裕度分析
[Gm,Pm,Wcg,Wcp] = margin(Go) 
Gm_dB = 20*log10(Gm)
margin(sys)

%% 求极点
s = tf('s')
G = 1/(s+1)
[num,den] = tfdata(G,'v') % 'v':以vector数组形式存储在[num,den]中
roots(den)

%% 从simulink中提取传递函数并分析伯德图
linsys1_tf = tf(linsys1)
linsys2_tf = tf(linsys2)
bode(linsys1_tf,linsys2_tf)
%...

%% state-space && transfer-function
linsys1_ss = ss(linsys2_tf) % tf -> ss
linsys1_tf = tf(linsys1_ss) % ss- > tf
```

