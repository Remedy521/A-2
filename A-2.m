clear;clc;close all;
global N ;%STAR elemnt number
N = 20;
global M ;%Antenna number of BS
M = 8;
NUM = 100; %信道次数
P_IM = -15:5:15;%IU发射功率dBm
PIM=10.^(P_IM/10);%转换成dBm
P_OM = -15:5:15;%IU发射功率dBm
POM=10.^(P_OM/10);%转换成dBm
Ka_Los = 5;
pno = sqrt(10.^(-105/10));
x_BS = 5i;
x_E = 0;
x_S = 50+10i;
x_IU = 50+15i;
x_OU = 50-15i;
q = 3;
for l = 1:NUM
% channel initialization
h_IS = sqrt(0.001*(abs(x_S-x_IU))^(-2.5))*1/sqrt(2)*(normrnd(0,1,N,1)+1i*normrnd(0,1,N,1));%IU到STAR信道
h_OS = sqrt(0.001*(abs(x_S-x_OU))^(-2.5))*1/sqrt(2)*(normrnd(0,1,N,1)+1i*normrnd(0,1,N,1));%OU到STAR信道
h_ES = sqrt(0.001*(abs(x_S-x_E))^(-2.5))/pno*1/sqrt(2)*(normrnd(0,1,N,1)+1i*normrnd(0,1,N,1));%STAR到E信道

G_NLoS = 1/sqrt(2)*(normrnd(0,1,N,M)+1i*normrnd(0,1,N,M));%STAR到BS地NLoS信道
G_LoS_t = 1;
G_LoS_r = 1;
d_r_lamda = 0.5;
ang_G_t = atan(0.1);
ang_G_r = pi-ang_G_t;
for i = 2:M
	G_LoS_t = [G_LoS_t,exp(1i*2*pi*(i-1)*d_r_lamda*sin(ang_G_t))];
end
for i = 2:N
	G_LoS_r = [G_LoS_r,exp(1i*2*pi*(i-1)*d_r_lamda*sin(ang_G_r))];
end
G_LoS = G_LoS_r'*G_LoS_t;

G = sqrt(0.001*(abs(x_BS-x_S))^(-2.2)*Ka_Los/(Ka_Los+1))*G_LoS+sqrt(0.001*(abs(x_BS-x_S))^(-2.2)*1/(Ka_Los+1))*G_NLoS;%STAR到BS信道
G = G/pno;

q_I = G'*diag(h_IS);
q_O = G'*diag(h_OS);

q_EI = h_ES'*diag(h_IS);
J_EI = q_EI'*q_EI;
q_EO = h_ES'*diag(h_OS);
J_EO = q_EO'*q_EO;

% power initialization
for np = 1:length(PIM)
PI = PIM(np);
PO = POM(np);

nd = 1;
Q_min_o = 100;
Q_min_j = 10;
while(1)
if(nd == 1)
w_r = 1/sqrt(2)*(normrnd(0,1,M,1)+1i*normrnd(0,1,M,1));
w_1 = w_r/norm(w_r);
W_1 = w_1'*w_1;

u_r_t = 1/sqrt(2)*(normrnd(0,1,N,1)+1i*normrnd(0,1,N,1));
for i = 1:N
	u_t1(i) = sqrt(0.5)*u_r_t(i)/abs(u_r_t(i));
end
U_t1 = u_t1'*u_t1;
u_r_r = 1/sqrt(2)*(normrnd(0,1,N,1)+1i*normrnd(0,1,N,1));
for i = 1:N
	u_r1(i) = sqrt(0.5)*u_r_r(i)/abs(u_r_r(i));
end
U_r1 = u_r1'*u_r1;

u_t_1 = u_t1.';
u_r_1 = u_r1.';

ow_1 = 1;
ns = 1;
lamda = 0;
end

I = eye(N);
I_2 = eye(M);

ro1 = 0.1;
ro2 = 0.1;
while(1)    
cvx_begin sdp
            variables gI_l gO_l gO_u Q z ow tau pn_1 pn_2 pn_3
            variable W(M,M) complex semidefinite
            variable U_t(N,N) complex semidefinite
            variable U_r(N,N) complex semidefinite
            variable u_t(N,1) complex
            variable u_r(N,1) complex
            maximize(tau-ro1*pn_1-ro2*pn_2-pn_3)
            subject to
                
                PI*gI_l >= PO*gO_u;
            
                U_t >= 0;
                U_r >= 0;
                for i = 1:N
                    U_t(i,i)+U_r(i,i) <= 1;
                end
                real(trace(U_t'*(I-u_t_1*u_t_1'))) <= pn_1;
                real(trace(U_r'*(I-u_r_1*u_r_1'))) <= pn_2;

                pn_1 >= 0;
                pn_2 >= 0;
                pn_3 >= 0;
                
                W >= 0;
                trace(W) == 1;
                real(trace(W'*(I_2-w_1'*w_1))) <= pn_3;
                
                2*real(trace((q_O'*W_1*q_O+U_r1)*(q_O'*W*q_O+U_r)'))-norm(q_O'*W_1*q_O+U_r1,'fro')^2 >= 4*gO_l+pow_pos(norm(q_O'*W*q_O-U_r,'fro'),2);
                4*gO_u+2*real(trace((q_O'*W_1*q_O-U_r1)*(q_O'*W*q_O-U_r)'))-norm(q_O'*W_1*q_O-U_r1,'fro')^2 >= pow_pos(norm(q_O'*W*q_O+U_r,'fro'),2);
                
                2*real(trace((q_I'*W_1*q_I+U_t1)*(q_I'*W*q_I+U_t)'))-norm(q_I'*W_1*q_I+U_t1,'fro')^2 >= 4*gI_l+pow_pos(norm(q_I'*W*q_I-U_t,'fro'),2);
                
                2*PI*gI_l >= pow_pos(Q*ow_1,2)+pow_pos((PO*gO_u+1)/ow_1,2)
               
                gO_l >= 0;
                Q >= 0;
                
                1+Q-lamda*(1+PI*real(trace(U_t*J_EI))) >= tau;
                1+PO*gO_l-lamda*(1+PO*real(trace(U_r*J_EO))) >= tau;
                
                
 cvx_end
 if cvx_optval==-Inf
	break;
 end 

    W_1 = W;
    U_t1 = U_t;
    U_r1 = U_r;
	ow_1 = sqrt((PO*gO_u+1)/Q);
    
    U_1 = U_t;
	[V_1 D_1] = eig(U_1);
  	D1= diag(D_1);
	[a1 b1] = max(D1);
	u_t_1 = V_1(:,b1);
    
    U_2 = U_r;
	[V_2 D_2] = eig(U_2);
  	D2= diag(D_2);
	[a2 b2] = max(D2);
	u_r_1 = V_2(:,b2);
    
    U_3 = W;
	[V_3 D_3] = eig(U_3);
  	D3= diag(D_3);
	[a3 b3] = max(D3);
	w_1 = V_3(:,b3);
    
    lamda1 = (1+Q)/(1+PI*real(trace(U_t*J_EI)));
    lamda2 = (1+PO*gO_l)/(1+PO*real(trace(U_r*J_EO)));
    lamda = min(lamda1,lamda2);
    
    R_s_conv_IU(ns,np,l) = log2(lamda);
    
    if(pn_1>=0.0001)
        ro1 = ro1*1.01;
    end
    if(pn_2>=0.0001)
        ro2 = ro2*1.01;
    end
    
	if(tau<=0.005)&&(max([pn_1,pn_2,pn_3])<=0.0001)&&(ns>2)
        u_tt = V_1*sqrt(D_1);
        u_t = u_tt(:,N);
        nt(np) = norm(u_t*u_t'-U_t);
        u_rr = V_2*sqrt(D_2);
        u_r = u_rr(:,N);
        nr(np) = norm(u_r*u_r'-U_r);
        ww = V_3*sqrt(D_3);
        w = ww(:,M);
        nw(np) = norm(w*w'-W);
        break;
    end
	ns = ns + 1;
end
if cvx_optval==-Inf
	break;
 end 
PI  = PIM(np);
A1 = real(trace(q_I'*W*q_I*U_t));
A2 = real(trace(U_t*J_EI));
B1 = real(trace(q_O'*W*q_O*U_r));
B2 = real(trace(U_r*J_EO));
s = 1;

R_sI_UW = log2((1+A1*PI/(B1*PO+1))/(1+A2*PI));
R_sO_UW = log2((1+B1*PO)/(1+B2*PO));

PO_tem = optimal_power(PI,A1,A2,B1,B2,s);
PO = min(PO_tem,POM(np));

R_sI_P = log2((1+A1*PI/(B1*PO+1))/(1+A2*PI));
R_sO_P = log2((1+B1*PO)/(1+B2*PO));



if mod(nd,2)==0
	Q_min_o = min(R_sI_P,R_sO_P);%n为偶数      
else
	Q_min_j = min(R_sI_P,R_sO_P);%n为奇数   
end
if(abs(Q_min_o-Q_min_j)<=0.001)
    
    R_c_s1(np,l) = max(Q_min_o,Q_min_j);
    
    %quantification
    [U_rq, u_rq] = Quantification(u_r,q,N);
    [U_tq, u_tq] = Quantification(u_t,q,N);
    A1 = real(trace(q_I'*W*q_I*U_tq));
    A2 = real(trace(U_tq*J_EI));
    B1 = real(trace(q_O'*W*q_O*U_rq));
    B2 = real(trace(U_rq*J_EO));
    s = 1;

    R_sI_q = log2((1+A1*PI/(B1*PO+1))/(1+A2*PI));
    R_sO_q = log2((1+B1*PO)/(1+B2*PO));
    
    R_q_s1(np,l) = min(R_sI_q,R_sO_q);
    clear R_sI_UW R_sO_UW R_sI_P R_sO_P R_sI_q R_sO_q;
    break;
end
nd = nd + 1;
end
if cvx_optval==-Inf
	break;
 end 

R_s_NOMA(np,l) = R_c_s1(np,l);

end
end
