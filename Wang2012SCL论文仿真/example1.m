%%  仿真说明
% 该仿真为王老师2012年论文的仿真验证与复现
%
% 题目："Observer design for discrete-time descriptor systems: An LMI approach"
%      “离散时间描述（广义）系统观测器设计：线性矩阵不等式方法”
% 期刊：《系统与控制快报》(Systems & Control Letters)
%

%%   仿真1：线性描述（广义）系统
% -------------------------------------------------------------------------
% Additional Function Needed: None
% Additional Toolbox Needed:  LMI Toolbox
% Additional Solver Needed:   None
% -------------------------------------------------------------------------
% Version:              1.0
% Author:               Weijie Ren
% Contact:              weijie.ren@outlook.com
% Initial modified:     Sep. 01, 2021
% Last modified:        
% -------------------------------------------------------------------------
% All rights reserved.
% Copyright (c) 2021, Weijie Ren. All rights reserved.
% Unauthorized copying of this file, via any medium, is strictly prohibited.
% -------------------------------------------------------------------------


%% 清空
clc;            % 清空命令行窗口
clear;          % 清空工作区
close all;      % 关闭所有图形视窗

%% 仿真参数
Nk = 200;       % 仿真总步数
Nks = 1:Nk;     % 序列生成（用于作图）
% Ts = ?;       % 采样时间，需要时可添加

%% 系统模型参数
E = [1 2 1; 0 2 1; 1 0 0];  % 描述矩阵
A = [0.153, 0.045, 0.069;   % 系统矩阵
     0.156, 0.252, 0.156;
     0.135, -0.171, -0.636];
B = [1; 1; 0.2];            % 输入矩阵
C = [1 0 0; 0 1 0];         % 输出矩阵

%% 系统维数
n = size(A,1);  % 状态向量x的维数
p = size(B,2);  % 输入向量u的维数
m = size(C,1);  % 输出向量y的维数

%% 系统等价变换矩阵计算
% 如下为文章中仿真1给出的T和N矩阵（错误！）
% 可取消注释测试其是否符合条件，同时注释掉下面正确的T和N
% T = [0 0 1; -0.4 -0.4 -0.2; 0.6 0.4 -0.8];
% N = [0 0; -0.2 1; 0.2 -2];

% 一组正确的T和N（手动计算结果）  % 注：该系统不太复杂，可以手工计算
T = [-1 1 1; -2 2 0; -2 3 -1];
N = [1 0; 2 1; 3 -2];
% 验证是否符合条件
rank(T)             % 检查方阵T是否满秩，即T的秩是否等于n
disp(T*E+N*C)       % 检查结果是否为单位阵In

%% 可观性/可检测性判断
% 常用的判别方法
rank(obsv(T*A,C))   % 若结果为n，则说明变换后的系统可观
                    % 即(TA,C)可观
                    % obsv()--计算可观判别式，参见现代控制理论书本

% 文中判别式1
rank([E; C])        % 是否等于n

% 文中判别式2
rank([1.5*E-A; C])  % 是否等于n，矩阵E前的系数z可取其他合适的值
                    % -以上三式可叠加使用-

%% 观测器增益求取
L = gain1(A,C,T);  % 获取相关信息可在命令行窗口输入 help gain1 或 doc gain1
% L = place((T*A)',C',[0.1; 0.2; 0.3])';
L
eig(T*A-L*C)    % 检验计算得到的增益L是否正确
                % 误差矩阵的所有特征值是否都位于单位圆内
                % 即 |特征值i|<1, i = 1,...,n

%% 预留变量内存
% 系统变量
x = zeros(n,Nk);    % 系统状态
u = zeros(p,Nk);    % 控制输入
y = zeros(m,Nk);    % 可测输出
z = zeros(1,Nk);    % 中间变量，求解此广义系统真实状态时所需
% 观测器变量
xhat = zeros(n,Nk); % 估计值
zeta = zeros(n,Nk); % 估计值中间变量

%% 初值
% 如下选取的初值和论文仿真1中几乎相同
x(:,1) = [1; 2; 0.4];               % 系统状态初值，仿真中可以任意给定
y(:,1) = C*x(:,1);                  % 初始输出值
xhat(:,1) = [0; 0; 0];              % 估计值初值，仿真中可以任意给定
zeta(:,1) = xhat(:,1) - N*y(:,1);   % 中间变量初值

%% 主循环
for k = 1:Nk
    % 系统控制输入（观测器设计时可任意给定输入）
    u(:,k) = 0.1*sin(0.1*k);
    
    % 系统真实状态（需要手动重新推导）
    if k>2
        x(2,k) = 61/126*z(k) - 23/189*x(1,k) - 100/567*u(:,k);
        x(3,k) = 46/183*x(1,k) + 4/61*x(2,k) + 200/549*u(:,k);
    end
    z(k+1) = A(2,:)*x(:,k) + u(:,k);
    x(1,k+1) = A(3,:)*x(:,k) + 0.2*u(:,k);
    % 系统输出
    y(:,k) = C*x(:,k);
    % 注：上式中部分系数因小数过多，故用分数表示，分数由科学计算器算得，其他表达亦可
    % 注：仿真时系统的真实状态值x可知，但对于实际系统，只能得到可测输出值y，内部的状态一般不可知
    
    % *特别注意*，若E不满秩（奇异）时，不可等式两边同时左乘E的伪逆
    % 错误示范：
    % x(:,k+1) = pinv(E)*A*x(:,k) + pinv(E)*B*u(:,k);
    
    % 观测器
    xhat(:,k) = zeta(:,k) + N*y(:,k);
    zeta(:,k+1) = T*A*xhat(:,k)+T*B*u(:,k)+L*(y(:,k)-C*xhat(:,k));
    % 中间变量选取 zeta = xhat - Ny
    % 文章中没有这个步骤，可以对照程序自己推导
end

%% 绘制图像
figure
plot(Nks,x(1,1:Nk),'b-')
hold on
plot(Nks,xhat(1,1:Nk),'r--')
xlabel('k/step')
ylabel('x_1(k)')
legend('x_1真实值','x_1估计值')

figure
plot(Nks,x(2,1:Nk),'b-')
hold on
plot(Nks,xhat(2,1:Nk),'r--')
xlabel('k/step')
ylabel('x_2(k)')
legend('x_2真实值','x_2估计值')

figure
plot(Nks,x(3,1:Nk),'b-')
hold on
plot(Nks,xhat(3,1:Nk),'r--')
xlabel('k/step')
ylabel('x_3(k)')
legend('x_3真实值','x_3估计值')
