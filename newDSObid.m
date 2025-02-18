clc
clear
warning('off');
%% 输入数据
% mpc=ukgds95;
mpc=bus2;
T=3;                                
% 基准值
Ub=mpc.bus(1,10);       
Sb=mpc.baseMVA;         
Zb=Ub^2/Sb;             
Ib=Ub/Zb;              
Yb=1/Zb;

n_bus=size(mpc.bus,1);          
n_bra=size(mpc.branch,1);       
n_gen=size(mpc.gen,1);         
n_dg=size(mpc.DGdata,1);


% 标幺化 
r_pu=mpc.branch(:,3);     %线路
x_pu=mpc.branch(:,4);     %线路
Gs_pu=mpc.bus(:,5);       %节点并联电导的标幺值
Bs_pu=mpc.bus(:,6);       %节点并联电纳的标幺值
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
c_dg2 = mpc.DGdata(:,4); %发电成本系数
c_dg1 = mpc.DGdata(:,5);%发电成本系数

c_curt2 = mpc.DGdata(:,7);%调整负荷成本系数
c_curt1 = mpc.DGdata(:,8);%调整负荷成本系数

c_str1 = mpc.DGdata(:,11);

P_strmax=mpc.DGdata(:,12);
cap_0=mpc.DGdata(:,13);
r_str=mpc.DGdata(:,14);
MSR = 0.003;%ramping capacity
R_pu=mpc.DGdata(:,9);%ramping up
R_pd=mpc.DGdata(:,10);%ramping down
E_exch = 0.5;%交换量

ARt=0.01;%FPP备用余量
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
%% 定义变量
% 节点支路数

% 定义变量 
xp_gen=sdpvar(n_gen,T);         % 节点注入有功
xq_gen=sdpvar(n_gen,T);         % 节点注入无功

xp_fbus=sdpvar(n_bra,T);        % 线路有功（指从fbus to tbus）
xq_fbus=sdpvar(n_bra,T);        % 线路无功（指从fbus to tbus）

xv_bus=sdpvar(n_bus,T);          % 节点电压的平方
xtheta_bus=sdpvar(n_bus,T);      % 节点电压的角度
xh_bra=sdpvar(n_bra,T);          % 支路电流的平方

xp_i_t=sdpvar(n_dg,T);           %分布式电源电能出力
xr_i_t=sdpvar(n_dg,T);           %分布式电源备用出力
xrdw_i_t=sdpvar(n_dg,T);  

xpcurt_i_t=sdpvar(n_dg,T);           %调整负荷电能出力
xrcurt_i_t=sdpvar(n_dg,T);           %调整负荷备用出力
pres_i_t=sdpvar(n_dg,T);

xpstr_i_t=sdpvar(n_dg,T);    %电池电能出力

E_t=sdpvar(1,T); %能量bid
R_t=sdpvar(1,T); %备用bid
Rdw_t=sdpvar(1,T);
P_P2P=sdpvar(n_dg,T);

d_ER_t=sdpvar(1,T); %FPP负荷
DGdata = mpc.DGdata;


%% 目标函数
%f=sum(r_pu.*xh_bra);
% f=sum(cdg*xp_dg+cg*xp_grid+xp_gen'*diag(c2)*xp_gen+c1'*xp_gen+c0);
%f=sum(xp_gen'*diag(c2)*xp_gen+c1'*xp_gen+c0);
% f=sum(sum(sum(c_dg2.*(xp_i_t+xr_i_t).^2+c_dg1.*(xp_i_t+xr_i_t)))+sum(sum(c_curt2.*(xpcurt_i_t+xrcurt_i_t).^2+c_curt1.*(xpcurt_i_t+xrcurt_i_t)))-sum(Pi_E_t.*E_t)-sum(Pi_R_t.*R_t));
% f=-(sum(Pi_E_t*E_t')+sum(Pi_R_t*R_t')-sum(sum(repmat(c_dg2,1,T).*(xp_i_t+xr_i_t).^2+repmat(c_dg1,1,T).*(xp_i_t+xr_i_t)))-sum(sum(repmat(c_str1,1,T).*(abs(xpstr_i_t)))) -sum(sum(repmat(c_curt2,1,T).*(xpcurt_i_t+xrcurt_i_t).^2+repmat(c_curt1,1,T).*(xpcurt_i_t+xrcurt_i_t))));
f=-(sum(Pi_E_t*E_t')+sum(Pi_R_t*(R_t'+Rdw_t')-sum(sum(repmat(c_dg2,1,T).*(xp_i_t).^2+repmat(c_dg1,1,T).*(xp_i_t)))-sum(sum(repmat(c_str1,1,T).*(abs(xpstr_i_t)))) -sum(sum(repmat(c_curt2,1,T).*(xpcurt_i_t).^2+repmat(c_curt1,1,T).*(xpcurt_i_t)))));

%% 约束条件 F=[F;...];
F=[];
% 节点*支路关联矩阵
Abus_bra=zeros(n_bus,n_bra);
for i=1:n_bra  
    Abus_bra(mpc.branch(i,1),i)=1;
    Abus_bra(mpc.branch(i,2),i)=-1;
end

% 节点与电源点的关联矩阵，存在电源点
Agen=zeros(n_bus,n_gen);
for i=1:n_gen
    Agen(mpc.gen(i,1),i)=1;
end

% 节点与流入功率的网损的关联矩阵
 Abus_bra2=Abus_bra;
 Abus_bra2(Abus_bra2==1)=0; % 确定所关联的节点（保留-1，即流入节点）
 Abus_bra2=-Abus_bra2;      % 确定网损为流出功率，正
 
 %节点与分布式电源关联矩阵
 Adg=zeros(n_bus,n_dg);
 for i=1:n_dg
     Adg(mpc.Prodata(i,1),i)=1;
 end

 
 %电压相角的约束
%  F=[F;Abus_bra'*xtheta_bus-(x_pu.*xp_fbus-r_pu.*xq_fbus)/(1*1)==zeros(n_bra,1);];

  for t=1:T
   for i=1:length(DGdata(:,1))
       F=[F;P_genmin(i,1)<=xp_i_t(i,t)<=P_genmax(i,1)];          %出力约束 
   end
  end
  
    for t=1:T
   for i=1:length(DGdata(:,1))
       F=[F;0<=pres_i_t(i,t)<=pres_max(i,t)];          %出力约束 
   end
  end
  
    for t=1:T
   for i=1:length(DGdata(:,1))
             F=[F;0<=xpcurt_i_t(i,t)<=P_curtmax(i,t)];          %负荷管理约束 
   end
    end
  
  for t=1:T
   for i=1:length(DGdata(:,1))
             F=[F;0<=xr_i_t(i,t)<=min(10*MSR,P_genmax(i,1)-xp_i_t(i,t))];
             F=[F;0<=xrdw_i_t(i,t)<=min(10*MSR,xp_i_t(i,t)-P_genmin(i,1))]; %备用最大约束              
   end
  end
  
  % for t=1:T 
  %     F=[F;sum(xp_i_t(:,t)+pres_i_t(:,t)-0.9*xpstr_i_t(:,t))-E_t(:,t)-d_ER_t(:,t)==0];   %供给平衡方程
  %      F=[F;d_ER_t(:,t)+sum(xpcurt_i_t(:,t))-D_t(:,t)==0]; 
  % end
  for t=1:T 
      F=[F;sum(xp_i_t(:,t)+pres_i_t(:,t)-P_P2P(:,t)-0.9*xpstr_i_t(:,t))-E_t(:,t)+sum(xpcurt_i_t(:,t))-D_t(:,t)==0];   %供给平衡方程
  end  

  for t=1:T 
      F=[F;0<=R_t(:,t)<=sum(xr_i_t(:,t)+pres_max(:,t)-pres_i_t(:,t))];   %备用平衡方程
      F=[F;0<=Rdw_t(:,t)<=sum(xrdw_i_t(:,t))];   %备用平衡方程
  end
  
   for t=2:T
   for i=1:length(DGdata(:,1))
             F=[F;xp_i_t(i,t)-xp_i_t(i,t-1)+xr_i_t(i,t)-xrdw_i_t(i,t-1)<=R_pu(i,:)];
             F=[F;xp_i_t(i,t-1)-xp_i_t(i,t)+xr_i_t(i,t-1)-xrdw_i_t(i,t)<=R_pd(i,:)]; %备用最大约束              
   end
  end

  
    for i=1:length(DGdata(:,1))
      F=[F;-cap_0<=sum(xpstr_i_t(i,:))<=(P_strmax-cap_0)];   %电池容量平衡方程
    end
    
  for t=1:T 
     for i=1:length(DGdata(:,1))
             F=[F;abs(xpstr_i_t(i,t))<=r_str(i,1)];          %电池充放电约束              
   end
  end
    
  
  for t=1:T 
      F=[F;abs(E_t)-E_exch<=0];
      F=[F;abs(E_t+R_t-Rdw_t)-E_exch<=0];   %交换电量约束
  end
  
  
%   for t=1:T 
%       F=[F;sum(P_genmax)-sum(xp_i_t(:,t)+xr_i_t(:,t))+sum(P_curtmax)-sum(xpcurt_i_t(:,t)+xrcurt_i_t(:,t))+E_exch-E_t(:,t)-R_t(:,t)-ARt>=0];   %备用余量方程
%   end
 for t=1:T
 F=[F;Abus_bra*xp_fbus(:,t)+(load_p(:,t)+Adg*(Load(:,t)-xpcurt_i_t(:,t)-xp_i_t(:,t)-pres_i_t(:,t)+0.9*xpstr_i_t(:,t)))-Agen*xp_gen(:,t)+Abus_bra2*(r_pu.*xh_bra(:,t))==zeros(n_bus,1);];     
% F=[F;Abus_bra*xp_fbus(:,t)+(load_p(:,t)+Adg*(Load(:,t)-xpcurt_i_t(:,t)-xrcurt_i_t(:,t)-xp_i_t(:,t)-xr_i_t(:,t)))+Gs_pu.*xv_bus(:,t)+Abus_bra2*(r_pu.*xh_bra(:,t))==zeros(n_bus,1);];     
 %无功平衡方程
 F=[F;Abus_bra*xq_fbus(:,t)+load_q(:,t)-Agen*xq_gen(:,t)+Abus_bra2*(x_pu.*xh_bra(:,t))==zeros(n_bus,1);];
 
 %支路电压降约束
 F=[F;Abus_bra'*xv_bus(:,t)-2*(r_pu.*xp_fbus(:,t)+x_pu.*xq_fbus(:,t))+(r_pu.^2+x_pu.^2).*xh_bra(:,t)==zeros(n_bra,1);];
 end
 %SOCP
 for t=1:T
 for i=1:n_bra
    F=[F;rcone([xp_fbus(i,t);xq_fbus(i,t)],0.5*xh_bra(i,t),xv_bus(mpc.branch(i,1),t))];  
%     F=[F;xp_fbus(i,1)^2+xq_fbus(i,1)^2==xh_bra(i,1)*xv_bus(mpc.branch(i,1),1)];
 end
 end
 

 %平衡节点电压约束和节点电压上下限边界约束
index_1 = find(mpc.bus(:,2)==3);
index_2 = find(mpc.bus(:,2)~=3);
for i=1:T
F = [F;xv_bus(index_1,i)==1,0.93^2<=xv_bus(index_2,i)<=1.07^2];
end
% F = [F;xv_bus(index_1,1)==1];      %平衡节点电压幅值

%% 求解
ops = sdpsettings('verbose',2,'solver', 'cplex','savesolveroutput',1);
result = solvesdp(F,f,ops);

%% 输出结果
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
        Current(i,t) = xh_bra(i,t)*Ib^2;%线路电流的平方
        Voltage(mpc.branch(i,1),t) = xv_bus(mpc.branch(i,1),t)*Ub^2;%线路电压
        Pfbus(i,t) = xp_fbus(i,t)*Sb;%线路有功功率
        Qfbus(i,t) = xq_fbus(i,t)*Sb;%线路无功功率
        error_socp(i,t) = Current(i,t)*Voltage(mpc.branch(i,1),t)-Pfbus(i,t)^2-Qfbus(i,t)^2;%各支路二阶锥规划的误差
    end
end

Max_error_socp = max(abs(error_socp));

f = double(f)
Max_error_socp = double(Max_error_socp);
E_t=double(E_t);bar(E_t*1000)
R_t=double(R_t);bar(R_t*1000)
