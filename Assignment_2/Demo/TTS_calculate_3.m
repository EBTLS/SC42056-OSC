function [TTS,y,x,z,m,s,u] = TTS_calculate_3(input)
% calculate the sum of TTS based on input u
u0=input(1);
u=input(2:end);
u_options=[15,20,25,30,35,40,45];
u0=u_options(u0);
u=u_options(u);
%%%%%%%%%%%%%%%%%%%%%%constant%%%%%%%%%%%%%%%%%%%%%%%%%%
c=60;
round=60;
l_veh=7;

E1=9;
E2=13;
E3=5;

%parameter of ud

N_ud=3;
v_ud=50;
l_ud=1000;

mu_udo1=1600/3600;
mu_udo2=1800/3600;
mu_udo3=1500/3600;

be_udo1=0.33;
be_udo2=0.34;
be_udo3=0.33;

c_ud=N_ud*l_ud/l_veh;

max_tau_ud=c_ud*l_veh/(N_ud*v_ud*c/3600*1000);

%parameter of o1d

N_o1d=3;
v_o1d=60;
l_o1d=1000;

mu_o1du=1600/3600;
mu_o1do3=1800/3600;
mu_o1do2=1500/3600;

be_o1du=0.33;
be_o1do2=0.34;
be_o1do3=0.33;

c_o1d=N_o1d*l_o1d/l_veh;

max_tau_o1d=c_o1d*l_veh/(N_o1d*v_o1d*60/3600*1000);


% time_varying "constant"
ae_ud=[(1800+10*E1)*ones(1,20),...
    (2100+10*E2)*ones(1,20),...
    (2300+10*E3)*ones(1,round-40)]/3600; %a_enter_ud

ae_o1d=(2000+10*E1)*ones(1,round)/3600; %a_enter_o1d

c_do1=zeros(1,round); 
for i=1:1:round
    if i<=20
        c_do1(i)=40+E1;
    elseif i<=35
        c_do1(i)=40+E1-2*(i-20);
    elseif i<=45
        c_do1(i)=10+E1;
    else
        c_do1(i)=10+E1+2*(i-45);
    end
end
c_do1;

c_do2=c_do1-E2*ones(1,round);

c_do3=[(30-E3)*ones(1,30),...
    (30+E3)*ones(1,round-30)];

c_du=[(40-E3)*ones(1,30),...
    (40+E3)*ones(1,round-30)];

C=[ae_ud;ae_o1d;c_du;c_do1;c_do2;c_do3];

C0=C(:,1);

%% states defination

%total states
x=zeros(10,1);  % states x
z=zeros(6,1);   % the time-varying "constant"
% u=30;   %input
m=zeros(6,1); % a^leave, the minimize function
s=zeros(4,1); %states i

tau=zeros(0,1);
lam=zeros(2,1); %first line for ud, second line for o1d

y=zeros(1,1);

%current states(initial)
temp_x=zeros(10,1);  % states x
temp_z=zeros(6,1);   % the time-varying "constant"
temp_u=u(1);   %input
temp_m=zeros(6,1); % a^leave, the minimize function
temp_s=zeros(4,1); %states i

temp_tau=zeros(2,1); %first line for ud, second line for o1d
temp_lam=zeros(2,1); %first line for ud, second line for o1d

temp_y=zeros(1,1);

%(k-1)states
prev_x=zeros(10,1);  % states x
prev_z=zeros(6,1);   % the time-varying "constant"
prev_u=30;   %input
prev_m=zeros(6,1); % a^leave, the minimize function
prev_s=zeros(4,1); % states i

prev_lam=zeros(2,1); %first line for ud, second line for o1d

prev_y=zeros(1,1);

%k-2 states
prev_prev_z=zeros(6,1); % the time-varying "constant"


%% initial variables(k=0)

%initial z
temp_z=C0;

%initial u
temp_u=u0;

%initial for lam & tau
temp_tau(1)=floor((c_ud-temp_s(3))*l_veh/ ...
    (N_ud*v_ud*c/3600*1000));
temp_lam(1)=mod((c_ud-temp_s(3))*l_veh/ ...
    (N_ud*v_ud/3600*1000),c);

temp_tau(2)=floor((c_o1d-temp_s(4))*l_veh/ ...
    (N_o1d*v_o1d*c/3600*1000));
temp_lam(2)=mod((c_o1d-temp_s(4))*l_veh/ ...
    (N_o1d*v_o1d/3600*1000),c);

%initial for aa_ud
if temp_tau(1)==0
    temp_x(3)=(c-temp_lam(1))/c*temp_z(1)+temp_lam(1)/c*prev_z(1);
    
elseif temp_tau(1)==1
    temp_x(3)=(c-temp_lam(1))/c*prev_z(1)+temp_lam(1)/c*prev_prev_z(1);
    
else
    "k-3 appear!";
    
end

%initial for aa_o1d
if temp_tau(2)==0
    temp_x(4)=(c-temp_lam(2))/c*temp_z(2)+temp_lam(2)/c*prev_z(2);

elseif temp_tau(1)==1
    temp_x(4)=(c-temp_lam(2))/c*prev_z(2)+temp_lam(2)/c*prev_prev_z(2);

else
    "k-3 appear!";
    
end

%initial al (mins)
temp_m(1)=min([mu_udo1*temp_u/c,...
    temp_x(5)/c+be_udo1*temp_x(3),...
    temp_z(4)/c/2]);
temp_m(2)=min([mu_udo2*temp_u/c,...
    temp_x(6)/c+be_udo2*temp_x(3),...
    temp_z(5)/c/2]);
temp_m(3)=min([mu_udo3*c/c,...
    temp_x(7)/c+be_udo3*temp_x(3),...
    temp_z(6)/c/2]);
temp_m(4)=min([mu_o1du*c/c,...
    temp_x(8)/c+be_o1du*temp_x(4),...
    temp_z(3)/c/2]);
temp_m(5)=min([mu_o1do2*(c-temp_u)/c,...
    temp_x(9)/c+be_o1do2*temp_x(4),...
    temp_z(5)/c/2]);
temp_m(6)=min([mu_o1do3*(c-temp_u)/c,...
    temp_x(10)/c+be_o1do3*temp_x(4),...
    temp_z(6)/c/2]);

temp_s(1)=temp_m(1)+temp_m(2)+temp_m(3);
temp_s(2)=temp_m(4)+temp_m(5)+temp_m(6);


x=temp_x;
y=temp_y;

prev_u=temp_u;
prev_x=temp_x;
prev_y=temp_y;
prev_z=temp_z;
prev_s=temp_s;
prev_mt=temp_m;
prev_tau=temp_tau;
prev_lam=temp_lam;



%% iteration
for k=1:1:60
    
    temp_z=C(:,k);
    temp_u=u(k);
    
    % iteration for n
    temp_x(1)=prev_x(1)+c*(prev_z(1)-prev_s(1));
    temp_x(2)=prev_x(2)+c*(prev_z(2)-prev_s(2));
    
    %iteration for y
    temp_y=c*(temp_x(1)+temp_x(2));
              
    %iteration for q
    temp_x(5)=prev_x(5)+(be_udo1*prev_x(3)-prev_m(1))*c;
%     prev_x(3)
%     prev_m(1)
%     be_udo1*prev_x(3)-prev_m(1)
    temp_x(6)=prev_x(6)+(be_udo2*prev_x(3)-prev_m(2))*c;
    temp_x(7)=prev_x(7)+(be_udo3*prev_x(3)-prev_m(3))*c;
    temp_x(8)=prev_x(8)+(be_o1du*prev_x(4)-prev_m(4))*c;
    temp_x(9)=prev_x(9)+(be_o1do2*prev_x(4)-prev_m(5))*c;
    temp_x(10)=prev_x(10)+(be_o1do3*prev_x(4)-prev_m(6))*c;
    
    %iteration for sum q
    temp_s(3)=temp_x(5)+temp_x(6)+temp_x(7);
    temp_s(4)=temp_x(8)+temp_x(9)+temp_x(10);
    
    %iteration for lam & tau 
    temp_tau(1)=floor((c_ud-temp_s(3))*l_veh/ ...
        (N_ud*v_ud*c/3600*1000));
    temp_lam(1)=mod((c_ud-temp_s(3))*l_veh/ ...
        (N_ud*v_ud/3600*1000),c);
    
    temp_tau(2)=floor((c_o1d-temp_s(4))*l_veh/ ...
        (N_o1d*v_o1d*c/3600*1000));
    temp_lam(2)=mod((c_o1d-temp_s(4))*l_veh/ ...
        (N_o1d*v_o1d/3600*1000),c);
    
    %iteration for aa_ud
    if temp_tau(1)==0
        temp_x(3)=(c-temp_lam(1))/c*temp_z(1)+temp_lam(1)/c*prev_z(1);       
    elseif temp_tau(1)==1
        temp_x(3)=(c-temp_lam(1))/c*prev_z(1)+temp_lam(1)/c*prev_prev_z(1);       
    else
        "k-3 appear!";       
    end
    
    %iteration for aa_o1d
    if temp_tau(2)==0
        temp_x(4)=(c-temp_lam(2))/c*temp_z(2)+temp_lam(2)/c*prev_z(2);   
    elseif temp_tau(1)==1
        temp_x(4)=(c-temp_lam(2))/c*prev_z(2)+temp_lam(2)/c*prev_prev_z(2);   
    else
        "k-3 appear!";       
    end
    
    %iteration al (mins)
    temp_m(1)=min([mu_udo1*temp_u/c,...
        temp_x(5)/c+be_udo1*temp_x(3),...
        temp_z(4)/c/2]);
    temp_m(2)=min([mu_udo2*temp_u/c,...
        temp_x(6)/c+be_udo2*temp_x(3),...
        temp_z(5)/c/2]);
    temp_m(3)=min([mu_udo3*c/c,...
        temp_x(7)/c+be_udo3*temp_x(3),...
        temp_z(6)/c/2]);
    temp_m(4)=min([mu_o1du*c/c,...
        temp_x(8)/c+be_o1du*temp_x(4),...
        temp_z(3)/c/2]);
    temp_m(5)=min([mu_o1do2*(c-temp_u)/c,...
        temp_x(9)/c+be_o1do2*temp_x(4),...
        temp_z(5)/c/2]);
    temp_m(6)=min([mu_o1do3*(c-temp_u)/c,...
        temp_x(10)/c+be_o1do3*temp_x(4),...
        temp_z(6)/c/2]);

    temp_s(1)=temp_m(1)+temp_m(2)+temp_m(3);
    temp_s(2)=temp_m(4)+temp_m(5)+temp_m(6);
    
    %contanate 
    x=[x,temp_x];
    z=[z,temp_z];
    m=[m,temp_m];
    s=[s,temp_s];    
    u=[u,temp_u];   
    y=[y,temp_y]; 
    
    %update, prepare for k+1
    prev_u=temp_u;
    prev_prev_z=prev_z;
    prev_y=temp_y;
    prev_x=temp_x;
    prev_m=temp_m;
    prev_z=temp_z;
    prev_s=temp_s;
    prev_tau=temp_tau;
    prev_lam=temp_lam;        
    
end
TTS=sum(y);

if nargout==1
    TTS;

    




end

