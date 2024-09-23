%% ����2�����������������壩ϵͳ
% �������1���ƣ������������������Ҫ��ע��

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

%% ���
clc;
clear;
close all;

%% �������
Nk = 200;
Nks = 1:Nk;
% Ts = ?;

%% ϵͳģ�Ͳ���
E = [1 0 0 -1; 0 1 0 0; 0 0 1 -1; 0 0 0 0];
A = [0.7 -1 0 0; 0 0.5 0 0; 0 1 0.3 0; 0 0 0 0];
B = [1 0; 0 1; 0 0; 0 0];
C = [0 0 1 0];

%% ϵͳά��
n = size(A,1);
p = size(B,2);
m = size(C,1);

%% ϵͳ�ȼ۱任�������
T = [1 0 -1 -2; 0 1 0 0; 0 0 0 -1; 0 0 -1 -1];
N = [1; 0; 1; 1];
% ��֤�Ƿ��������
rank(T)             % T�Ƿ����ȣ�T���Ƿ����n
disp(T*E+N*C)       % ����Ƿ�Ϊ��λ��In

%% �ɹ���/�ɼ�����ж�
% ���õ��б𷽷�
rank(obsv(T*A,C))   % �����Ϊn����˵���任���ϵͳ�ɹ�
                    % ����������Tȡ��ֵ���÷�������ʾϵͳ����ȫ�ɹ�

% �����б�ʽ1
rank([E; C])        % �Ƿ����n

% �����б�ʽ2
rank([1.5*E-A; C])  % �Ƿ����n������Eǰ��ϵ��z��ȡ�������ʵ�ֵ

%% �۲���������ȡ
gamma = 0.05;   % ����ϣ�ģ�Lipschitz������
L = gain2(A,C,T,gamma);  % ��ȡ�����Ϣ���������д������� help gain2 �� doc gain2
L
eig(T*A-L*C)
% *��Ҫ���벻Ҫ������ϣ�ĳ���gamma��Ϊ������LMI�����㣬gamma����Ϊ����*

%% Ԥ�������ڴ�
% ϵͳ����
x = zeros(n,Nk);
u = zeros(p,Nk);
y = zeros(m,Nk);
% �۲�������
xhat = zeros(n,Nk);
zeta = zeros(n,Nk);

%% ��ֵ
x(:,1) = [-1; 2; 1; 0.5];
y(:,1) = C*x(:,1);
xhat(:,1) = [0; 0; 0; 0];
zeta(:,1) = xhat(:,1) - N*y(:,1);

%% ��ѭ��
for k = 1:Nk
    % ��������
    u(:,k) = [0.1*sin(0.1*k);
              0];
    
    % ϵͳ��ʵ״̬
    x(4,k+1) = 0.5*sin(0.1*k+1);  % ϵͳ״̬4Ϊ�����ɱ�����
    x(1,k+1) = A(1,:)*x(:,k) + u(1,k) + x(4,k+1);
    x(2,k+1) = A(2,:)*x(:,k) + u(2,k) + 0.05*sin(x(3,k));
    x(3,k+1) = A(3,:)*x(:,k) + x(4,k+1);
    % ϵͳ���
    y(:,k) = C*x(:,k);
    
    % �۲���
    xhat(:,k) = zeta(:,k) + N*y(:,k);
    zeta(:,k+1) = T*A*xhat(:,k)+T*B*u(:,k)+L*(y(:,k)-C*xhat(:,k))...
                  +T*[0; 0.05*sin(xhat(3,k)); 0; 0];
end

%% ����ͼ��
figure
plot(Nks,x(1,1:Nk),'b-')
hold on
plot(Nks,xhat(1,1:Nk),'r--')
xlabel('k/step')
ylabel('x_1(k)')
legend('x_1��ʵֵ','x_1����ֵ')

figure
plot(Nks,x(2,1:Nk),'b-')
hold on
plot(Nks,xhat(2,1:Nk),'r--')
xlabel('k/step')
ylabel('x_2(k)')
legend('x_2��ʵֵ','x_2����ֵ')

figure
plot(Nks,x(3,1:Nk),'b-')
hold on
plot(Nks,xhat(3,1:Nk),'r--')
xlabel('k/step')
ylabel('x_3(k)')
legend('x_3��ʵֵ','x_3����ֵ')

figure
plot(Nks,x(4,1:Nk),'b-')
hold on
plot(Nks,xhat(4,1:Nk),'r--')
xlabel('k/step')
ylabel('x_4(k)')
legend('x_4��ʵֵ','x_4����ֵ')
