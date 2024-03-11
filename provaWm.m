s=tf('s');
k=ureal('k',0.09,'Range',[0.089 0.4]);
tk=k.NominalValue;


b=ureal('b',0.019,'Range',[0.0036 0.04]);
tb=b.NominalValue;

gn=(10*tb*s + 10*tk)/(s^2 * (s^2 + 11*tb*s + 11*tk));
gn
x1=(0.076*s +4)/(s^2 * (s^2 + 0.0836*s + 4.4));
x0=(0.19*s+0.9)/(s^2*(s^2 + 0.209*s +0.99));
x10= x1-x0;
x10o=inv(x0)*x10;


%g0=(0.19*s+0.9)/(s^2*(s^2 + 0.209*s +0.99));
%g0=(10*b.NominalValue*s +10*k.NominalValue)/(s^2*(s^2 + 11*b.NominalValue*s + 11*k.NominalValue));
sd=(s^2 +0.44*s)/(s^2 +0.44*s +0.16);
td=0.16/(s^2 +0.44*s +0.16);
W1=makeweight(100,0.2,0.5);
W3=makeweight(0.3,1,100);
W2=(s+5)/(s+20);
Wm=(1.8*x10o); % quanto vale?
 delta=ultidyn('delta',[1 1]);
 
g11=gn*(1 + delta*Wm);
g11.InputName='u';
g11.OutputName='yg';
W1.InputName='y';
W1.OutputName='z1';
W2.InputName='u';
W2.OutputName='z2';
W3.InputName='yg';
W3.OutputName='z3';
Sumy=sumblk('y=w-yg');
P=connect(g11,W1,W2,W3,Sumy,{'w','u'},{'z1','z2','z3','y'});


 [K,clperf,info]=musyn(P,1,1);
 clperf
K=minreal(K);
bode(K);
g11=minreal(g11);
looptransfer=loopsens(g11,K); 
%  Ti=looptransfer.Ti;
%  [PeakNorm,freq]=norm(Ti.NominalValue,'inf'); % torna picco massimo
% % %peaknorm deve essere minore di 1!! altrimenti non è stab robusto
% % % verifica stab robusta con valori singolari strutturati
%  PeakNorm
To=looptransfer.To;
step(To);
eig(To)
 omega=logspace(-1,2,500);
 opt=robOptions('Display','on');
 [stabmarg,wcu,info]=robstab(To,omega,opt);
 stabmarg
P=minreal(P);
 CLP=lft(P,K);
[PERFMAG,WCU_PERF,INFO]=robgain(CLP,1.5,omega,opt);
PERFMAG


  omega=logspace(-4,1,1000); 
figure(1)
sigma(inv(W1),'r-',sd,'b--',omega)
axis([10^(-4) 10^1 -100 10])
grid
title(' confronto W1 e Sd','FontSize',24)
legend('1/W1','sd')
ax=gca;
ax.FontSize=20;
  omega2=logspace(-2,4,1000); % se non esce stabile, ovvero in robsta

figure(2)
sigma(inv(W3),'r-',td,'b--',omega2)
axis([10^(-2) 10^4 -160 20])
grid
title(' confronto W3 e Td','FontSize',24)
legend('1/W3','td')
ax=gca;
ax.FontSize=20;

omega3=logspace(-3,3,1000);
figure(3)
sigma(inv(W3),'r-',td,'b--',To,'y-',omega3)
axis([10^(-3) 10^3 -250 50])
grid
title(' confronto 1/W3,Td e To','FontSize',24)
legend('1/W3','td','To')
ax=gca;
ax.FontSize=20;
omega4=logspace(-4,3,1000);
figure(4)
sigma(inv(W1),'r-',sd,'b--',looptransfer.So,'y-',omega4)
axis([10^(-4) 10^3 -90 10])
grid
title(' confronto 1/W1,Sd e So','FontSize',24)
legend('1/W1','Sd','So')
ax=gca;
ax.FontSize=20;


