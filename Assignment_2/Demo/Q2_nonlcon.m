function [c,ceq] = Q2_nonlcon(x)

%the nonlinear constraints about Q2

%% constant
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

ae_o1d=(2100+10*E1)*ones(1,round)/3600; %a_enter_o1d

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

%% iteration

for k=1:1:round
    if (k==1)
    
    else
        (x(25+(k-1)*31)-(mu_udo1*x(31+(k-1)*31)))...
            *(x(25+(k-1)*31)-(x(5+(k-1)*31)/c+be_udo1*x(3+(k-1)*31)))...
            *(x(25+(k-1)*31)-(x(14+(k-1)*31)/c/2));
        (x(26+(k-1)*31)-(mu_udo2*x(31+(k-1)*31)))...
            *(x(26+(k-1)*31)-(x(6+(k-1)*31)/c+be_udo2*x(3+(k-1)*31)))...
            *(x(26+(k-1)*31)-(x(15+(k-1)*31)/c/2));
        (x(27+(k-1)*31)-(mu_udo3))...
            *(x(27+(k-1)*31)-(x(7+(k-1)*31)/c+be_udo3*x(3+(k-1)*31)))...
            *(x(27+(k-1)*31)-(x(16+(k-1)*31)/c/2));
        (x(28+(k-1)*31)-(mu_o1du))...
            *(x(28+(k-1)*31)-(x(8+(k-1)*31)/c+be_o1du*x(4+(k-1)*31)))...
            *(x(28+(k-1)*31)-(x(13+(k-1)*31)/c/2));
        (x(29+(k-1)*31)-(mu_o1do2*(c-x(31+(k-1)*31))/c))...
            *(x(29+(k-1)*31)-(x(9+(k-1)*31)/c+be_o1do2*x(4+(k-1)*31)))...
            *(x(29+(k-1)*31)-(x(15+(k-1)*31)/c/2));
        (x(30+(k-1)*31)-(mu_o1do3*(c-x(31+(k-1)*31))/c))...
            *(x(29+(k-1)*31)-(x(10+(k-1)*31)/c+be_o1do3*x(4+(k-1)*31)))...
            *(x(29+(k-1)*31)-(x(16+(k-1)*31)/c/2));
       
    end
    

end

ceq=x(1)-(x(2)/2)*(x(2)>=x(3))-(x(3)/2)*(x(2)<x(3))

end

