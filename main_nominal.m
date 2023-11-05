clear all
close all

global J k1 delta 

J=[1 -1 0; -1 5 0; 0 0 7];
k1=2;
delta=3*k1;%(0,4k1)

% q_0=[1;0;0;0];
% o_0=[0;0;0];
% h_0=1;

q_0=[-1;0;0;0];
o_0=[10;10;10];
h_0=1;

xi_0=[q_0;o_0;h_0];

TSPAN = [0 50];
JSPAN = [0 100];
rule  = 1;

options = odeset('AbsTol',1e-6,'RelTol',1e-3,'InitialStep',eps);
    
[t,j,xi] = HyEQsolver(@F,@G,@C,@D,xi_0,TSPAN,JSPAN,rule,options);

save data_nominal 

%q(1)-h vs 0
%o vs 0
q_vec=xi(:,1:4);
o_vec=xi(:,5:7);
h_vec=xi(:,8);
q1_vec=xi(:,1);
figure(1)
plot(t,abs(q1_vec-h_vec))
set(gca, 'YScale', 'log')
figure(2)
plot(t,sum(o_vec.^2,2).^(0.5))
set(gca, 'YScale', 'log')
figure(3)
plot(t,h_vec)
set(gca, 'YScale', 'linear')


% %check if in A_xi={V(z)\leq\eigmax(P)*c^2} eventually====>plot V(z) and the line \eigmax(P)*c^2
% V_vec=zeros(length(t),1);
% for i=1:1:length(t)
%     V_vec(i)=fcn_V(transpose(xi(i,:)));
% end
% c=1/theta*(min(eigs(K)))^(-1)*alpha*norm(delta);
% plot(t,V_vec,t,max(eigs(P))*c^2.*ones(length(t),1))
% set(gca, 'YScale', 'log')

save("data_nominal","q_vec","o_vec","h_vec","-append")

function dxi = F(xi)
    global J k1 delta
    q=xi(1:4);
    o=xi(5:7);
    h=xi(8);
    
    u=fcn_kappa(xi);
    
    dq=1/2.*fcn_E(q)*o;
    do=J\(fcn_S(J*o)*o+u);
    dh=0;
    dxi = [dq;do;dh];
end

function next_xi = G(xi)
    global J k1 delta
    q=xi(1:4);
    o=xi(5:7);
    h=xi(8);
    
    next_xi = [q;o;-h];
end

function out = C(xi)
    global J k1 delta 
    
    out = (fcn_mu(xi)<=delta);
end

function out = D(xi)
    global J k1 delta
    
    out = (fcn_mu(xi)>=delta);
end

function out=fcn_kappa(in)
    global J k1 delta
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    out=-k1.*h.*q(2:4)-o;%theta(o)
end

function out=fcn_mu(in)
    global J k1 delta
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
    global J k1 delta
    q=in(1:4);
    o=in(5:7);
    h=in(8);
    
    out=2*k1*(1-h*q(1))+1/2.*o'*J*o;
end