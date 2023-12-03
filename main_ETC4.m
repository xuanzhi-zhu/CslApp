clear all
close all

rng('shuffle')

global J k1 delta sigma c_lb_prime_prime chi_ub

J=[1 -1 0; -1 5 0; 0 0 7];
k1=2;
delta=0.01*k1;%(0,4k1)
sigma=-1;%(0,1)
chi_ub=1;

% c_lb=0.1;%>0

% c_lb_prime=2*k1+0.5*delta+0.5*max(eigs(J))*c_lb^2;

% c_lb_prime_prime=c_lb_prime;

% q_0=[1;0;0;0];
% o_0=[0;0;0];
% h_0=1;
% hu_0=[0;0;0];

% %can have complete discrete sol
% q_0=[-1;0;0;0];
% o_0=[0;0;0];
% h_0=1;
% hu_0=[0;0;0];

%
% aux=[2;9;-7;2];
% q_0=aux./norm(aux);
% q_0=-[1;0;0;0];
q_0=[-0.1;0.7;0.5;0.5];
% o_0=[10;-10;10];
o_0=[0;0;0];
h_0=1;
hu_0=[-1;1;-1];
counter_0=0;
chi_0=1;

xi_0=[q_0;o_0;h_0;hu_0;counter_0;chi_0];

TSPAN = [0 1];
JSPAN = [0 200];
rule  = 1;%1 for jumps
% ,'InitialStep',eps
options = odeset('AbsTol',1e-6,'RelTol',1e-3,'InitialStep',eps);

[t,j,xi] = HyEQsolver(@F,@G,@C,@D,xi_0,TSPAN,JSPAN,rule,options);

save data_ETC4


q_vec=xi(:,1:4);
o_vec=xi(:,5:7);
h_vec=xi(:,8);
hu_vec=xi(:,9:11);
q1_vec=xi(:,1);
% figure(1)
% plot(t,abs(q1_vec-h_vec))
% set(gca, 'YScale', 'log')
% figure(2)
% plot(t,sum(o_vec.^2,2).^(0.5))
% set(gca, 'YScale', 'log')
% figure(3)
% plot(t,h_vec)
% set(gca, 'YScale', 'linear')
% figure(4)
% plot(t,hu_vec)
% set(gca, 'YScale', 'linear')

%check ETC
C_etc=zeros(length(t),1);
C_syn=zeros(length(t),1);
C_vec=zeros(length(t),1);
for i=1:1:length(t)
    aux=transpose(xi(i,:));
    gradV=fcn_gradV(aux);
    f=F(aux);
    f=f(1:11);
    W=fcn_W(aux);
%     o_aux=transpose(o_vec(i,:));
    
    chi_aux=xi(i,13);
    
    C_etc(i)=(gradV'*f<=sigma*W)&(chi_aux<=chi_ub);
    C_syn(i)=(fcn_mu(aux)<=delta);
    C_vec(i)=C_etc(i)&C_syn(i);
end


%check ETC
D_etc=zeros(length(t),1);
D_syn=zeros(length(t),1);
D_vec=zeros(length(t),1);
for i=1:1:length(t)
    aux=transpose(xi(i,:));
    gradV=fcn_gradV(aux);
    f=F(aux);
    f=f(1:11);
    W=fcn_W(aux);
%     o_aux=transpose(o_vec(i,:));
    chi_aux=xi(i,13);    

    D_etc(i)=(gradV'*f>=sigma*W)&(chi_aux>=chi_ub)&(chi_aux<=chi_ub+eps);
    D_syn(i)=(fcn_mu(aux)>=delta);
    D_vec(i)=D_etc(i)|D_syn(i);
end
t_etc=t(D_etc==1);
t_etc=t_etc(1:end);

t_inter_etc=t_etc-[0;t_etc(1:end-1)];%inter-event time
t_inter_etc=t_inter_etc(1:end);

figure(1)
plot(t_etc,t_inter_etc,'x');
set(gca, 'YScale', 'log')


figure(2)
%check if x eventually in A s.t. V(x)<=c_lb_prime
V_vec=zeros(length(t),1);
for i=1:1:length(t)
    chi_aux=xi(13);
    V_vec(i)=exp(1*chi_aux)*abs(fcn_V(transpose(xi(i,:))));
end
% plot(t,V_vec,t,c_lb_prime_prime.*ones(length(t),1))
plot(t,V_vec)

set(gca, 'YScale', 'log')



figure(3)

C_=zeros(length(t),1);
D_=zeros(length(t),1);
for i=1:1:length(t)
    aux=transpose(xi(i,:));
    C_(i)=C(aux);
    D_(i)=D(aux);
end


% save("data_ETC4","V_vec","c_lb_prime_prime","t_etc","t_inter_etc","q_vec","o_vec","h_vec","hu_vec","D_etc","-append")

function dxi = F(xi)
    global J k1 delta sigma c_lb_prime_prime
    q=xi(1:4);
    o=xi(5:7);
    h=xi(8);
    hu=xi(9:11);
    
    u=hu;
    
    dq=1/2.*fcn_E(q)*o;
    do=J\(fcn_S(J*o)*o+u);
    dh=0;
    dhu=[0;0;0];
    
    dxi = [dq;do;dh;dhu;0;1];
end

function next_xi = G(xi)
    global J k1 delta sigma chi_ub
    q=xi(1:4);
    o=xi(5:7);
    h=xi(8);
    hu=xi(9:11);
    counter=xi(12);
    chi=xi(13);
    
    gradV=fcn_gradV(xi);
    f=F(xi);
    f=f(1:11);
    W=fcn_W(xi);
    
    
    
    D_syner=(fcn_mu(xi)>=delta);
    next_syner=[q;o;-h;hu;counter;chi];
    
    D_etc=(gradV'*f>=sigma*W)&(chi>=chi_ub)&(chi<=chi_ub+eps);
    next_etc=[q;o;h;fcn_kappa([q;o;h]);counter+1;0];
    
%     if counter==0
%         next_xi=next_etc;
%     elseif counter==1
%         next_xi=next_syner;
%     elseif counter==2
%         next_xi=next_etc;
%     else
        I2=round(1+rand(1,1));%1 or 2
        next_syner_etc=[next_syner next_etc];

        if (D_syner==1)&&(D_etc==0)
            next_xi = next_syner;
        elseif (D_syner==0)&&(D_etc==1)
            next_xi = next_etc;
        else
            next_xi = next_syner_etc(:,2);%to check if ETC can issue complete discrete sols, change I2 to 2
        end
%     end
end

function out = C(xi)
    global J k1 delta sigma chi_ub
%     o=xi(5:7);
    
    gradV=fcn_gradV(xi);
    f=F(xi);
    f=f(1:11);
    W=fcn_W(xi);
    counter=xi(12);
    chi=xi(13);
    
    C_syner=(fcn_mu(xi)<=delta);
    C_etc=(gradV'*f<=sigma*W)&(chi<=chi_ub+eps);%%%%%%
    
%     if counter<=3
%         out=0;
%     else
        out = C_syner & C_etc ;
%     end
end

function out = D(xi)
    global J k1 delta sigma chi_ub
%     o=xi(5:7);
    
    gradV=fcn_gradV(xi);
    f=F(xi);
    f=f(1:11);
    W=fcn_W(xi);
    counter=xi(12);
    chi=xi(13);
    
    aux1=fcn_mu(xi);
    aux2=gradV'*f;
    aux3=sigma*W;
    
    D_syner=(fcn_mu(xi)>=delta);
    D_etc=(gradV'*f>=sigma*W)&(chi>=chi_ub)&(chi<=chi_ub+eps);
    
    out = (D_syner | D_etc)&(chi<=chi_ub+eps);
end

function out=fcn_kappa(in)
    global J k1 delta sigma
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    out=-k1.*h.*q(2:4)-o;%theta(o)
end

function out=fcn_mu(in)
    global J k1 delta sigma
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    x1=in;
    x2=[q;o;-h];
    
    out=fcn_V(x1)-min(fcn_V(x1),fcn_V(x2));
end

function out=fcn_R(q)
    out=quat2rotm(q');
end

function out=fcn_E(q)
    n=q(1);
    e=q(2:4);
    out=[-e';n*eye(3)+fcn_S(e)];
end

function out=fcn_S(x)
    out=[0 -x(3) x(2);...
        x(3) 0 -x(1);...
        -x(2) x(1) 0];
end

function out=fcn_V(in)
    global J k1 delta sigma
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    out=2*k1*(1-h*q(1))+1/2.*o'*J*o;
end

function out=fcn_gradV(in)
    global J k1 delta sigma
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    out=[-2.*k1.*h.*[1;0;0;0];J*o;0;zeros(3,1)];
end

function out=fcn_W(in)
    global J k1 delta sigma
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    out=-o'*o;%theta(o)
end