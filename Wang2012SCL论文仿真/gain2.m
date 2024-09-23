function L = gain2(A,C,T,gamma)
% GAIN2 函数用于求取-非线性-观测器增益矩阵L
%   该求解方法基于线性矩阵不等式(LMI)技术
%     函数调用格式：L = GAIN2(A,C,T,gamma)
%           输入部分
%               A--系统矩阵
%               C--输出矩阵
%               T--系统等价变换矩阵
%               gamma--李普希茨(Lipschitz)常数
%           输出部分
%               L--非线性观测器增益矩阵
%
%   Designed by WJ Ren, 1 September, 2021


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

%% 参数获取
n = size(A,1);
m = size(C,1);
gamma1 = gamma*max(eig(T));

%% 初始化LMI系统
setlmis([])

%% 定义决策变量
P = lmivar(1,[n 1]);
W = lmivar(2,[n m]);
eta = lmivar(1,[1 0]);  % 标量

%% 添加LMI内因子项
% #1 (P>0, P正定)
lmiterm([-1 1 1 P],1,1)
% #2 (eta>0)
lmiterm([-2 1 1 eta],1,1)
% #3
lmiterm([-3 1 1 P],-1,1)
lmiterm([-3 1 1 eta],1,gamma1^2)
lmiterm([-3 2 1 P],1,T*A)
lmiterm([-3 2 1 W],-1,C)
lmiterm([-3 2 2 P],1,1)
lmiterm([-3 2 2 eta],-1,1)
lmiterm([-3 3 1 P],1,T*A)
lmiterm([-3 3 1 W],-1,C)
lmiterm([-3 3 2 0],0)     % 内因子为零时可以不描述，描述了也没错
lmiterm([-3 3 3 P],-1,1)

%% 获取LMI系统描述
lmisys = getlmis;

%% LMI解算
[tmin,xfeas] = feasp(lmisys);

%% 变量操作（向量决策变量→矩阵变量）
P = dec2mat(lmisys,xfeas,P);
W = dec2mat(lmisys,xfeas,W);

L = P\W;  % 解得观测器增益L

end