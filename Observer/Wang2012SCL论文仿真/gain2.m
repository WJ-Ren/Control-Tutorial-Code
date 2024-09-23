function L = gain2(A,C,T,gamma)
% GAIN2 ����������ȡ-������-�۲����������L
%   ����ⷽ���������Ծ��󲻵�ʽ(LMI)����
%     �������ø�ʽ��L = GAIN2(A,C,T,gamma)
%           ���벿��
%               A--ϵͳ����
%               C--�������
%               T--ϵͳ�ȼ۱任����
%               gamma--����ϣ��(Lipschitz)����
%           �������
%               L--�����Թ۲����������
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

%% ������ȡ
n = size(A,1);
m = size(C,1);
gamma1 = gamma*max(eig(T));

%% ��ʼ��LMIϵͳ
setlmis([])

%% ������߱���
P = lmivar(1,[n 1]);
W = lmivar(2,[n m]);
eta = lmivar(1,[1 0]);  % ����

%% ���LMI��������
% #1 (P>0, P����)
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
lmiterm([-3 3 2 0],0)     % ������Ϊ��ʱ���Բ�������������Ҳû��
lmiterm([-3 3 3 P],-1,1)

%% ��ȡLMIϵͳ����
lmisys = getlmis;

%% LMI����
[tmin,xfeas] = feasp(lmisys);

%% �����������������߱��������������
P = dec2mat(lmisys,xfeas,P);
W = dec2mat(lmisys,xfeas,W);

L = P\W;  % ��ù۲�������L

end