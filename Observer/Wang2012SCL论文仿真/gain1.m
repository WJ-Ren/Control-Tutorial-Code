function L = gain1(A,C,T)
% GAIN1 ����������ȡ-����-�۲����������L
%   ����ⷽ���������Ծ��󲻵�ʽ(LMI)����
%     �������ø�ʽ��L = GAIN1(A,C,T)
%           ���벿��
%               A--ϵͳ����
%               C--�������
%               T--ϵͳ�ȼ۱任����
%           �������
%               L--���Թ۲����������
%
%   Designed by WJ Ren, 1 September, 2021

% �������Ծ��󲻵�ʽ��matlab��⣬�ɲο��ٷ�LMI Toolbox�ĵ���ֱ��help��doc+���ָ��


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

%% ��ʼ��LMIϵͳ
setlmis([])

%% ������߱���
P = lmivar(1,[n 1]);
W = lmivar(2,[n m]);

%% ���LMI��������
% #1 (P>0, P����)
lmiterm([-1 1 1 P],1,1)
% #2
lmiterm([2 1 1 P],-1,1)
lmiterm([2 2 1 P],1,T*A)
lmiterm([2 2 1 W],-1,C)
lmiterm([2 2 2 P],-1,1)

%% ��ȡLMIϵͳ����
lmisys = getlmis;

%% LMI����
[tmin,xfeas] = feasp(lmisys);

%% �����������������߱��������������
P = dec2mat(lmisys,xfeas,P);
W = dec2mat(lmisys,xfeas,W);

L = P\W;  % ��ù۲�������L

end