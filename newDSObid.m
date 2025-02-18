clc
clear
warning('off');
%% ��������
% mpc=ukgds95;
mpc=bus2;
T=3;                                
% ��׼ֵ
Ub=mpc.bus(1,10);       
Sb=mpc.baseMVA;         
Zb=Ub^2/Sb;             
Ib=Ub/Zb;              
Yb=1/Zb;

n_bus=size(mpc.bus,1);          
n_bra=size(mpc.branch,1);       
n_gen=size(mpc.gen,1);         
n_dg=size(mpc.DGdata,1);


% ���ۻ� 
r_pu=mpc.branch(:,3);     %��·
x_pu=mpc.branch(:,4);     %��·
Gs_pu=mpc.bus(:,5);       %�ڵ㲢���絼�ı���ֵ
Bs_pu=mpc.bus(:,6);       %�ڵ㲢�����ɵı���ֵ
load_p=zeros(n_bus,T);
% load_p=repmat(mpc.bus(:,3)/Sb,1,T);
load_q=repmat(mpc.bus(:,4)/Sb,1,T);
D_t=zeros(1,T);
Load=zeros(n_dg,T);
for t=1:T
for i=1:length(mpc.Prodata(:,1))
Load(i,t)=mpc.Prodata(i,t+1);
end
D_t(:,t)=sum(Load(:,t));
end

pres_max=zeros(n_dg,T);
for t=1:T
for i=1:length(mpc.RESdata(:,1))
pres_max(i,t)=mpc.RESdata(i,t+1);
end
end

P_curtmax=zeros(n_dg,T);
for t=1:T
for i=1:length(mpc.Pcurtmaxdata(:,1))
P_curtmax(i,t)=mpc.Pcurtmaxdata(i,t+1)/4;
end
end

P_genmin=mpc.DGdata(:,2);
P_genmax=mpc.DGdata(:,3);
c_dg2 = mpc.DGdata(:,4); %����ɱ�ϵ��
c_dg1 = mpc.DGdata(:,5);%����ɱ�ϵ��

c_curt2 = mpc.DGdata(:,7);%�������ɳɱ�ϵ��
c_curt1 = mpc.DGdata(:,8);%�������ɳɱ�ϵ��

c_str1 = mpc.DGdata(:,11);

P_strmax=mpc.DGdata(:,12);
cap_0=mpc.DGdata(:,13);
r_str=mpc.DGdata(:,14);
MSR = 0.003;%ramping capacity
R_pu=mpc.DGdata(:,9);%ramping up
R_pd=mpc.DGdata(:,10);%ramping down
E_exch = 0.5;%������

ARt=0.01;%FPP��������
A=30;
Bid=zeros(A,3);

 Pi_E_t=zeros(1,T);
 for t=1:T
Pi_E_t(1,t)=mpc.PriceE(1,t)*1000;
 end
 Pi_R_t=zeros(1,T);
 for t=1:T
Pi_R_t(1,t)=mpc.PriceR(1,t)*1000;
 end

% for a=1:A
% % Pi_R_t=4000+400*a;
% Pi_R_t=0;
% % Pi_E_t=10000;
% Pi_E_t=4000+333*a;
%% �������
% �ڵ�֧·��

% ������� 
xp_gen=sdpvar(n_gen,T);         % �ڵ�ע���й�
xq_gen=sdpvar(n_gen,T);         % �ڵ�ע���޹�

xp_fbus=sdpvar(n_bra,T);        % ��·�й���ָ��fbus to tbus��
xq_fbus=sdpvar(n_bra,T);        % ��·�޹���ָ��fbus to tbus��

xv_bus=sdpvar(n_bus,T);          % �ڵ��ѹ��ƽ��
xtheta_bus=sdpvar(n_bus,T);      % �ڵ��ѹ�ĽǶ�
xh_bra=sdpvar(n_bra,T);          % ֧·������ƽ��

xp_i_t=sdpvar(n_dg,T);           %�ֲ�ʽ��Դ���ܳ���
xr_i_t=sdpvar(n_dg,T);           %�ֲ�ʽ��Դ���ó���
xrdw_i_t=sdpvar(n_dg,T);  

xpcurt_i_t=sdpvar(n_dg,T);           %�������ɵ��ܳ���
xrcurt_i_t=sdpvar(n_dg,T);           %�������ɱ��ó���
pres_i_t=sdpvar(n_dg,T);

xpstr_i_t=sdpvar(n_dg,T);    %��ص��ܳ���

E_t=sdpvar(1,T); %����bid
R_t=sdpvar(1,T); %����bid
Rdw_t=sdpvar(1,T);
P_P2P=sdpvar(n_dg,T);

d_ER_t=sdpvar(1,T); %FPP����
DGdata = mpc.DGdata;


%% Ŀ�꺯��
%f=sum(r_pu.*xh_bra);
% f=sum(cdg*xp_dg+cg*xp_grid+xp_gen'*diag(c2)*xp_gen+c1'*xp_gen+c0);
%f=sum(xp_gen'*diag(c2)*xp_gen+c1'*xp_gen+c0);
% f=sum(sum(sum(c_dg2.*(xp_i_t+xr_i_t).^2+c_dg1.*(xp_i_t+xr_i_t)))+sum(sum(c_curt2.*(xpcurt_i_t+xrcurt_i_t).^2+c_curt1.*(xpcurt_i_t+xrcurt_i_t)))-sum(Pi_E_t.*E_t)-sum(Pi_R_t.*R_t));
% f=-(sum(Pi_E_t*E_t')+sum(Pi_R_t*R_t')-sum(sum(repmat(c_dg2,1,T).*(xp_i_t+xr_i_t).^2+repmat(c_dg1,1,T).*(xp_i_t+xr_i_t)))-sum(sum(repmat(c_str1,1,T).*(abs(xpstr_i_t)))) -sum(sum(repmat(c_curt2,1,T).*(xpcurt_i_t+xrcurt_i_t).^2+repmat(c_curt1,1,T).*(xpcurt_i_t+xrcurt_i_t))));
f=-(sum(Pi_E_t*E_t')+sum(Pi_R_t*(R_t'+Rdw_t')-sum(sum(repmat(c_dg2,1,T).*(xp_i_t).^2+repmat(c_dg1,1,T).*(xp_i_t)))-sum(sum(repmat(c_str1,1,T).*(abs(xpstr_i_t)))) -sum(sum(repmat(c_curt2,1,T).*(xpcurt_i_t).^2+repmat(c_curt1,1,T).*(xpcurt_i_t)))));

%% Լ������ F=[F;...];
F=[];
% �ڵ�*֧·��������
Abus_bra=zeros(n_bus,n_bra);
for i=1:n_bra  
    Abus_bra(mpc.branch(i,1),i)=1;
    Abus_bra(mpc.branch(i,2),i)=-1;
end

% �ڵ����Դ��Ĺ������󣬴��ڵ�Դ��
Agen=zeros(n_bus,n_gen);
for i=1:n_gen
    Agen(mpc.gen(i,1),i)=1;
end

% �ڵ������빦�ʵ�����Ĺ�������
 Abus_bra2=Abus_bra;
 Abus_bra2(Abus_bra2==1)=0; % ȷ���������Ľڵ㣨����-1��������ڵ㣩
 Abus_bra2=-Abus_bra2;      % ȷ������Ϊ�������ʣ���
 
 %�ڵ���ֲ�ʽ��Դ��������
 Adg=zeros(n_bus,n_dg);
 for i=1:n_dg
     Adg(mpc.Prodata(i,1),i)=1;
 end

 
 %��ѹ��ǵ�Լ��
%  F=[F;Abus_bra'*xtheta_bus-(x_pu.*xp_fbus-r_pu.*xq_fbus)/(1*1)==zeros(n_bra,1);];

  for t=1:T
   for i=1:length(DGdata(:,1))
       F=[F;P_genmin(i,1)<=xp_i_t(i,t)<=P_genmax(i,1)];          %����Լ�� 
   end
  end
  
    for t=1:T
   for i=1:length(DGdata(:,1))
       F=[F;0<=pres_i_t(i,t)<=pres_max(i,t)];          %����Լ�� 
   end
  end
  
    for t=1:T
   for i=1:length(DGdata(:,1))
             F=[F;0<=xpcurt_i_t(i,t)<=P_curtmax(i,t)];          %���ɹ���Լ�� 
   end
    end
  
  for t=1:T
   for i=1:length(DGdata(:,1))
             F=[F;0<=xr_i_t(i,t)<=min(10*MSR,P_genmax(i,1)-xp_i_t(i,t))];
             F=[F;0<=xrdw_i_t(i,t)<=min(10*MSR,xp_i_t(i,t)-P_genmin(i,1))]; %�������Լ��              
   end
  end
  
  % for t=1:T 
  %     F=[F;sum(xp_i_t(:,t)+pres_i_t(:,t)-0.9*xpstr_i_t(:,t))-E_t(:,t)-d_ER_t(:,t)==0];   %����ƽ�ⷽ��
  %      F=[F;d_ER_t(:,t)+sum(xpcurt_i_t(:,t))-D_t(:,t)==0]; 
  % end
  for t=1:T 
      F=[F;sum(xp_i_t(:,t)+pres_i_t(:,t)-P_P2P(:,t)-0.9*xpstr_i_t(:,t))-E_t(:,t)+sum(xpcurt_i_t(:,t))-D_t(:,t)==0];   %����ƽ�ⷽ��
  end  

  for t=1:T 
      F=[F;0<=R_t(:,t)<=sum(xr_i_t(:,t)+pres_max(:,t)-pres_i_t(:,t))];   %����ƽ�ⷽ��
      F=[F;0<=Rdw_t(:,t)<=sum(xrdw_i_t(:,t))];   %����ƽ�ⷽ��
  end
  
   for t=2:T
   for i=1:length(DGdata(:,1))
             F=[F;xp_i_t(i,t)-xp_i_t(i,t-1)+xr_i_t(i,t)-xrdw_i_t(i,t-1)<=R_pu(i,:)];
             F=[F;xp_i_t(i,t-1)-xp_i_t(i,t)+xr_i_t(i,t-1)-xrdw_i_t(i,t)<=R_pd(i,:)]; %�������Լ��              
   end
  end

  
    for i=1:length(DGdata(:,1))
      F=[F;-cap_0<=sum(xpstr_i_t(i,:))<=(P_strmax-cap_0)];   %�������ƽ�ⷽ��
    end
    
  for t=1:T 
     for i=1:length(DGdata(:,1))
             F=[F;abs(xpstr_i_t(i,t))<=r_str(i,1)];          %��س�ŵ�Լ��              
   end
  end
    
  
  for t=1:T 
      F=[F;abs(E_t)-E_exch<=0];
      F=[F;abs(E_t+R_t-Rdw_t)-E_exch<=0];   %��������Լ��
  end
  
  
%   for t=1:T 
%       F=[F;sum(P_genmax)-sum(xp_i_t(:,t)+xr_i_t(:,t))+sum(P_curtmax)-sum(xpcurt_i_t(:,t)+xrcurt_i_t(:,t))+E_exch-E_t(:,t)-R_t(:,t)-ARt>=0];   %������������
%   end
 for t=1:T
 F=[F;Abus_bra*xp_fbus(:,t)+(load_p(:,t)+Adg*(Load(:,t)-xpcurt_i_t(:,t)-xp_i_t(:,t)-pres_i_t(:,t)+0.9*xpstr_i_t(:,t)))-Agen*xp_gen(:,t)+Abus_bra2*(r_pu.*xh_bra(:,t))==zeros(n_bus,1);];     
% F=[F;Abus_bra*xp_fbus(:,t)+(load_p(:,t)+Adg*(Load(:,t)-xpcurt_i_t(:,t)-xrcurt_i_t(:,t)-xp_i_t(:,t)-xr_i_t(:,t)))+Gs_pu.*xv_bus(:,t)+Abus_bra2*(r_pu.*xh_bra(:,t))==zeros(n_bus,1);];     
 %�޹�ƽ�ⷽ��
 F=[F;Abus_bra*xq_fbus(:,t)+load_q(:,t)-Agen*xq_gen(:,t)+Abus_bra2*(x_pu.*xh_bra(:,t))==zeros(n_bus,1);];
 
 %֧·��ѹ��Լ��
 F=[F;Abus_bra'*xv_bus(:,t)-2*(r_pu.*xp_fbus(:,t)+x_pu.*xq_fbus(:,t))+(r_pu.^2+x_pu.^2).*xh_bra(:,t)==zeros(n_bra,1);];
 end
 %SOCP
 for t=1:T
 for i=1:n_bra
    F=[F;rcone([xp_fbus(i,t);xq_fbus(i,t)],0.5*xh_bra(i,t),xv_bus(mpc.branch(i,1),t))];  
%     F=[F;xp_fbus(i,1)^2+xq_fbus(i,1)^2==xh_bra(i,1)*xv_bus(mpc.branch(i,1),1)];
 end
 end
 

 %ƽ��ڵ��ѹԼ���ͽڵ��ѹ�����ޱ߽�Լ��
index_1 = find(mpc.bus(:,2)==3);
index_2 = find(mpc.bus(:,2)~=3);
for i=1:T
F = [F;xv_bus(index_1,i)==1,0.93^2<=xv_bus(index_2,i)<=1.07^2];
end
% F = [F;xv_bus(index_1,1)==1];      %ƽ��ڵ��ѹ��ֵ

%% ���
ops = sdpsettings('verbose',2,'solver', 'cplex','savesolveroutput',1);
result = solvesdp(F,f,ops);

%% ������
xp_gen  = double(xp_gen);
xq_gen  = double(xq_gen);
xp_fbus = double(xp_fbus);
xq_fbus = double(xq_fbus);
xv_bus  = double(xv_bus);
xtheta_bus=double(xtheta_bus);
xtheta_bus=xtheta_bus*180/pi;
xp_i_t=double(xp_i_t);
xr_i_t=double(xr_i_t);
xpcurt_i_t=double(xpcurt_i_t);
xpstr_i_t=double(xpstr_i_t);

Vbus = sqrt(xv_bus);
Ibra = sqrt(xh_bra);

for t=1:T
    for i=1:n_bra
        Current(i,t) = xh_bra(i,t)*Ib^2;%��·������ƽ��
        Voltage(mpc.branch(i,1),t) = xv_bus(mpc.branch(i,1),t)*Ub^2;%��·��ѹ
        Pfbus(i,t) = xp_fbus(i,t)*Sb;%��·�й�����
        Qfbus(i,t) = xq_fbus(i,t)*Sb;%��·�޹�����
        error_socp(i,t) = Current(i,t)*Voltage(mpc.branch(i,1),t)-Pfbus(i,t)^2-Qfbus(i,t)^2;%��֧·����׶�滮�����
    end
end

Max_error_socp = max(abs(error_socp));

f = double(f)
Max_error_socp = double(Max_error_socp);
E_t=double(E_t);bar(E_t*1000)
R_t=double(R_t);bar(R_t*1000)
