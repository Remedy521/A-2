function [PO] = optimal_power(P1,A1,A2,B1,B2,s)
%OPTIMAL_POWER �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
PO = (2*A2*B1*P1*s - A1*B2*P1*s + B1*s^2 - B2*s^2 - sqrt((-2*A2*B1*P1*s + A1*B2*P1*s - B1*s^2 + B2*s^2)^2 - 4*(-A2*B1^2*P1 - B1^2*s + B1*B2*s)*(A1*P1*s^2 - A2*P1*s^2)))/(2*(-A2*B1^2*P1 - B1^2*s + B1*B2*s));
end

