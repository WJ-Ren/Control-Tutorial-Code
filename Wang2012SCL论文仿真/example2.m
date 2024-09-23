%% 仿真2：非线性描述（广义）系统
% 因与仿真1类似，本仿真仅给出部分重要的注释

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
clc;
clear;
close all;

%% 仿真参数
Nk = 200;
Nks = 1:Nk;
% Ts = ?;

%% 系统模型参数
E = [1 0 0 -1; 0 1 0 0; 0 0 1 -1; 0 0 0 0];
A = [0.7 -1 0 0; 0 0.5 0 0; 0 1 0.3 0; 0 0 0 0];
B = [1 0; 0 1; 0 0; 0 0];
C = [0 0 1 0];

%% 系统维数
n = size(A,1);
p = size(B,2);
m = size(C,1);

%% 系统等价变换矩阵计算
T = [1 0 -1 -2; 0 1 0 0; 0 0 0 -1; 0 0 -1 -1];
N = [1; 0; 1; 1];
% 验证是否符合条件
rank(T)             % T是否满秩，T秩是否等于n
disp(T*E+N*C)       % 结果是否为单位阵In

%% 可观性/可检测性判断
% 常用的判别方法
rank(obsv(T*A,C))   % 若结果为n，则说明变换后的系统可观
                    % 本例中无论T取何值，该方法皆显示系统不完全可观

% 文中判别式1
rank([E; C])        % 是否等于n

% 文中判别式2
rank([1.5*E-A; C])  % 是否等于n，矩阵E前的系数z可取其他合适的值

%% 观测器增益求取
gamma = 0.05;   % 李普希茨（Lipschitz）常数
L = gain2(A,C,T,gamma);  % 获取相关信息可在命令行窗口输入 help gain2 或 doc gain2
L
eig(T*A-L*C)
% *重要！请不要将李普希茨常数gamma视为变量用LMI来计算，gamma请人为给定*

%% 预留变量内存
% 系统变量
x = zeros(n,Nk);
u = zeros(p,Nk);
y = zeros(m,Nk);
% 观测器变量
xhat = zeros(n,Nk);
zeta = zeros(n,Nk);

%% 初值
x(:,1) = [-1; 2; 1; 0.5];
y(:,1) = C*x(:,1);
xhat(:,1) = [0; 0; 0; 0];
zeta(:,1) = xhat(:,1) - N*y(:,1);

%% 主循环
for k = 1:Nk
    % 控制输入
    u(:,k) = [0.1*sin(0.1*k);
              0];
    
    % 系统真实状态
    x(4,k+1) = 0.5*sin(0.1*k+1);  % 系统状态4为“自由变量”
    x(1,k+1) = A(1,:)*x(:,k) + u(1,k) + x(4,k+1);
    x(2,k+1) = A(2,:)*x(:,k) + u(2,k) + 0.05*sin(x(3,k));
    x(3,k+1) = A(3,:)*x(:,k) + x(4,k+1);
    % 系统输出
    y(:,k) = C*x(:,k);
    
    % 观测器
    xhat(:,k) = zeta(:,k) + N*y(:,k);
    zeta(:,k+1) = T*A*xhat(:,k)+T*B*u(:,k)+L*(y(:,k)-C*xhat(:,k))...
                  +T*[0; 0.05*sin(xhat(3,k)); 0; 0];
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

figure
plot(Nks,x(4,1:Nk),'b-')
hold on
plot(Nks,xhat(4,1:Nk),'r--')
xlabel('k/step')
ylabel('x_4(k)')
legend('x_4真实值','x_4估计值')
