% The function used are defined in the final part of this .m file

%% global constant
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
    C0=C(:,1);

    u_starting_point_15 = 15*ones(1,61);
    u_starting_point_45 = 45*ones(1,61);

%% Task 3
%%%%%%%%%%%%%%%%%%%%%%Method 3: Matlab Optimization Toolbox

    fun=@TTS_calculate;
    options_sqp=optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',5000,'PlotFcn','optimplotfval');
    options_interior_point=optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',5000,'PlotFcn','optimplotfval');
    options_newton=optimoptions('fminunc','HessUpdate','dfp','MaxFunctionEvaluations',5000,'PlotFcn','optimplotfval');

    [x_1,fval_1,exitflag_1,output_1]=fmincon(fun,u_starting_point_15,[],[],[],[],15*ones(1,61),45*ones(1,61),[],options_sqp);
    [x_2,fval_2,exitflag_2,output_2]=fmincon(fun,u_starting_point_45,[],[],[],[],15*ones(1,61),45*ones(1,61),[],options_sqp);
    [x_3,fval_3,exitflag_3,output_3]=fmincon(fun,u_starting_point_15,[],[],[],[],15*ones(1,61),45*ones(1,61),[],options_interior_point);
    [x_4,fval_4,exitflag_4,output_4]=fmincon(fun,u_starting_point_45,[],[],[],[],15*ones(1,61),45*ones(1,61),[],options_interior_point);
    [x_5,fval_5,exitflag_5,output_5]=fminunc(fun,u_starting_point_15,options_newton);
    [x_6,fval_6,exitflag_6,output_6]=fminunc(fun,u_starting_point_45,options_newton);
    
    [TTS_sqp_1,y_sqp_1,x_sqp_1,z_sqp_1,m_sqp_1,s_sqp_1,u_sqp_1] = TTS_calculate(x_1);    
    [TTS_sqp_2,y_sqp_2,x_sqp_2,z_sqp_2,m_sqp_2,s_sqp_2,u_sqp_2] = TTS_calculate(x_2); 
    
    [TTS_inter_1,y_inter_1,x_inter_1,z_inter_1,m_inter_1,s_inter_1,u_inter_1] = TTS_calculate(x_3);    
    [TTS_inter_2,y_inter_2,x_inter_2,z_inter_2,m_inter_2,s_inter_2,u_inter_2] = TTS_calculate(x_4); 
    
    [TTS_newton_1,y_newton_1,x_newton_1,z_newton_1,m_newton_1,s_newton_1,u_newton_1] = TTS_calculate(x_5);    
    [TTS_newton_2,y_newton_2,x_newton_2,z_newton_2,m_newton_2,s_newton_2,u_newton_2] = TTS_calculate(x_6); 
    
    disp("SQP methods from 15s:");
    disp("optimal TTS:")
    disp(fval_1);
    disp("flag:")
    disp(exitflag_1)
    disp("iteration times")
    disp(output_1.iterations)
    
    disp("SQP methods from 45s:");
    disp("optimal TTS:")
    disp(fval_2);
    disp("flag:")
    disp(exitflag_2)
    disp("iteration times")
    disp(output_2.iterations)
    
    disp("Interior-point methods from 15s:");
    disp("optimal TTS:")
    disp(fval_3);
    disp("flag:")
    disp(exitflag_3)
    disp("iteration times")
    disp(output_3.iterations)
    
    disp("Interior-point methods from 45s:");
    disp("optimal TTS:")
    disp(fval_4);
    disp("flag:")
    disp(exitflag_4)
    disp("iteration times")
    disp(output_4.iterations)
    
    disp("Quasi-Newton methods from 15s:");
    disp("optimal TTS:")
    disp(fval_5);
    disp("flag:")
    disp(exitflag_5)
    disp("iteration times")
    disp(output_5.iterations)
    
    disp("Quasi-Newton methods from 45s:");
    disp("optimal TTS:")
    disp(fval_6);
    disp("flag:")
    disp(exitflag_6)
    disp("iteration times")
    disp(output_6.iterations)
    

%%%%%%%%%%%%%%%%%%%%%%Method 2: Simuannealing Optimizaiton Method

    options_simul=optimoptions('simulannealbnd','FunctionTolerance',1,'PlotFcn','saplotbestf','PlotInterval',100);

    %initial point

    [u_sim_1,fval_sim_1,exitflag_sim_1,output_sim_1]=simulannealbnd(fun,u_starting_point_15,15*ones(1,61),45*ones(1,61),options_simul);
    [u_sim_2,fval_sim_2,exitflag_sim_2,output_sim_2]=simulannealbnd(fun,u_starting_point_45,15*ones(1,61),45*ones(1,61),options_simul);

    disp("simul methods from 15s:");
    disp("optimal TTS:")
    disp(fval_sim_1);
    disp("flag:")
    disp(exitflag_sim_1)
    disp("iteration times")
    disp(output_sim_1.iterations)

    disp("simul methods from 45s:");
    disp("optimal TTS");
    disp(fval_sim_2);
    disp("flag:")
    disp(exitflag_sim_2)
    disp("iteration times")
    disp(output_sim_2.iterations)

    [TTS_sim_1,y_sim_1,x_sim_1,z_sim_1,m_sim_1,s_sim_1,u_sim_1] = TTS_calculate(u_sim_1);
    [TTS_sim_2,y_sim_2,x_sim_2,z_sim_2,m_sim_2,s_sim_2,u_sim_2] = TTS_calculate(u_sim_2);  

%Besides, we also tries the result with perpendicular direction search
%%%%%%%%%%%%%%%%%%%%%%Method 3: Perpendicular direction Search+Golden Section+Simuannealing Optimizaiton Method

%initial value of input
%for starting point 15s

    u_starting_point_15 = 15*ones(1,61);
    u_starting_point_45 = 45*ones(1,61);
    iter = 20;               
    count = 0;
    golden_iter = 3;


TTS_previous_2_2_1 = TTS_calculate(u_starting_point_15);
[u_gold_1, TTS_current_2_2_1,TTS_history_2_2_1] =...
    goldensection(u_starting_point_15,TTS_previous_2_2_1,golden_iter);
disp("optimal solution (15s) after perpendicular search with golden section steps:");
disp(TTS_current_2_2_1);
disp("iteration times")
disp(length(TTS_history_2_2_1))

%for starting point 45s

TTS_previous_2_2_2 = TTS_calculate(u_starting_point_45);
[u_gold_2, TTS_current_2_2_2,TTS_history_2_2_2] =...
    goldensection(u_starting_point_45,TTS_previous_2_2_2,golden_iter);
disp("optimal solution (45s) after perpendicular search with golden section steps:");
disp(TTS_current_2_2_2);
disp("iteration times")
disp(length(TTS_history_2_2_2))

[TTS_gold_1,y_gold_1,x_gold_1,z_gold_1,m_gold_1,s_gold_1,u_gold_1] = TTS_calculate(u_gold_1);    
[TTS_gold_2,y_gold_2,x_gold_2,z_gold_2,m_gold_2,s_gold_2,u_gold_2] = TTS_calculate(u_gold_2);    



%% Task 4
    u_starting_point_30=30*ones(1,61);
    [TTS_ref,y_ref,x_ref,z_ref,m_ref,s_ref,u_ref] = TTS_calculate(u_starting_point_30);
    
    figure('Name',"Task 4 traffic situation on control -2 each lane， u(t)=30s")

        subplot(2,3,1)
        bar([0:1:round],x_ref(5,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d,o1)")

        subplot(2,3,2)
        bar([0:1:round],x_ref(6,:),'b');
        hold on;
        grid on;
        title("number of vehicles,lane(u,d,o2)")

        subplot(2,3,3)
        bar([0:1:round],x_ref(7,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d,o3)")

        subplot(2,3,4)
        bar([0:1:round],x_ref(8,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d,u)")
        
        subplot(2,3,5)
        bar([0:1:round],x_ref(9,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d,o2)")

        subplot(2,3,6)
        bar([0:1:round],x_ref(10,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d,o3)")
        
    figure('Name',"Task 4 traffic situation SQP -2 each lane， u(t)=45s")

        subplot(2,3,1)
        bar([0:1:round],x_sqp_1(5,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d,o1)")

        subplot(2,3,2)
        bar([0:1:round],x_sqp_1(6,:),'b');
        hold on;
        grid on;
        title("number of vehicles,lane(u,d,o2)")

        subplot(2,3,3)
        bar([0:1:round],x_sqp_1(7,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d,o3)")

        subplot(2,3,4)
        bar([0:1:round],x_sqp_1(8,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d,u)")
        
        subplot(2,3,5)
        bar([0:1:round],x_sqp_1(9,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d,o2)")

        subplot(2,3,6)
        bar([0:1:round],x_sqp_1(10,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d,o3)")
    
    
    figure('Name',"Task 4 traffic situation on control， u(t)=30s")

        subplot(2,2,1)
        bar([0:1:round],x_ref(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_ref(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,lane(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],x_ref(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_ref(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")
        
    figure('Name','Task 4 traffic situation SQP,15s')
        subplot(2,2,1)
        bar([0:1:round],x_sqp_1(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_sqp_1(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_sqp_1(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_sqp_1(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")
        
    figure('Name','Task 4 traffic situation SQP,45s')
        subplot(2,2,1)
        bar([0:1:round],x_sqp_2(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_sqp_2(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_sqp_2(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_sqp_2(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")

    figure('Name','Task 4 traffic situation simulannealing,15s')
        subplot(2,2,1)
        bar([0:1:round],x_sim_1(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_sim_1(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_sim_1(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_sim_1(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")
    
    figure('Name','Task 4 traffic situation simulannealing,45s')
        subplot(2,2,1)
        bar([0:1:round],x_sim_2(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_sim_2(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_sim_2(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_sim_2(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")

    figure('Name','Task 4 traffic situation perpendicualr,15s')
        subplot(2,2,1)
        bar([0:1:round],x_gold_1(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_gold_1(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_gold_1(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_gold_1(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")
    
    figure('Name','Task 4 traffic situation perpendicualr,45s')
        subplot(2,2,1)
        bar([0:1:round],x_gold_2(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_gold_2(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_gold_2(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_gold_2(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")
        
     
    figure('Name','Task 4 input Matlab Tools')
        subplot(2,3,1)
        bar([0:1:round],x_1);
        hold on;
        grid on;
        title("sqp:15")
        
        subplot(2,3,4)
        bar([0:1:round],x_2);
        hold on;
        grid on;
        title("sqp:45")
 
        subplot(2,3,2)
        bar([0:1:round],x_3);
        hold on;
        grid on;
        title("interior:15")
        
        subplot(2,3,5)
        bar([0:1:round],x_4);
        hold on;
        grid on;
        title("interior:45")
        
        subplot(2,3,3)
        bar([0:1:round],x_5);
        hold on;
        grid on;
        title("Newton 15")
        
        subplot(2,3,6)
        bar([0:1:round],x_1);
        hold on;
        grid on;
        title("Newton 45")
        

    figure('Name','Task 4 perpendicular+golden_section process')

        subplot(2,2,1)
        plot(TTS_history_2_2_1);
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_history_2_2_1)),'r');
        title("(15s):process")
        legend('optimized','u(t)=30s')

        subplot(2,2,3)
        bar(u_gold_1);
        grid on;
        title("(15s):u(t)")

        subplot(2,2,2)
        plot(TTS_history_2_2_2);
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_history_2_2_2)),'r');
        title(" (45s):process")
        legend('optimized','u(t)=30s')

        subplot(2,2,4)
        bar(u_gold_2);
        grid on;
        title(" (15s):u(t)")

%% Task 5
    %%%genetic optimization method

        %discrete green light time set
        u_options=[15,20,25,30,35,40,45];

        %algorithm settings
        generation_round=300;
        member_number=100;
        %the exchange_position of parents
        exchange_position=30;
        %mutation probabilitTTS
        mutation_probability=0.2;
        %parents number
        parents_number=20;

        %%random initial generation

        %initial generation 0
        u_generation_5_1_1=randi(length(u_options),[1,round]);
        %u_generation=4*ones(member_number,round);
        u0_generation_5_1_1=randi(length(u_options),[member_number,1]);
        for i =2:1:member_number
            u_generation_5_1_1=[u_generation_5_1_1;randi(length(u_options),[1,round])];
        end

        [u_final_5_1_1,TTS_optimal_history_5_1_1]=...
            genetic_optimization(u0_generation_5_1_1,u_generation_5_1_1,u_options,...
                round,generation_round,member_number,mutation_probability,exchange_position,parents_number);
        TTS_optimal_5_1_1=TTS_optimal_history_5_1_1(end);
        u_final_5_1_1=(u_final_5_1_1+2)*5;
        disp("TTS of genetic+initiaal points based on random initial points 1");
        disp(TTS_optimal_5_1_1);

        %%%%%%%%%initial generation 2

        %initial generation 0
        u_generation_5_1_2=randi(length(u_options),[1,round]);
        u0_generation_5_1_2=randi(length(u_options),[member_number,1]);
        for i =2:1:member_number
            u_generation_5_1_2=[u_generation_5_1_2;randi(length(u_options),[1,round])];
        end

        for i=1:1:length(u_options)
                u_generation_5_1_2(i,:)=i*ones(1,round);
                u0_generation_5_1_2(i)=i; 
        end
        [u_final_5_1_2,TTS_optimal_history_5_1_2]=...
            genetic_optimization(u0_generation_5_1_2,u_generation_5_1_2,u_options,...
                round,generation_round,member_number,mutation_probability,exchange_position,parents_number);
        TTS_optimal_5_1_2=TTS_optimal_history_5_1_2(end);
        u_final_5_1_2=(u_final_5_1_2+2)*5;
        disp("TTS of genetic+initiaal points based on random initial points 2");
        disp(TTS_optimal_5_1_2);

        %%%%%%%%%initial generation based on 30s
        %initial generation 
        u_generation_5_1_3=4*ones(member_number,round);
        u0_generation_5_1_3=4*ones(member_number,1);

        for i=1:1:member_number
            for j=1:1:round
            poss=rand(1);
            if (poss<=0.2)
                u_generation_5_1_3(i,j)=randi(length(u_options),1);
            end
            end
        end

        [u_final_5_1_3,TTS_optimal_history_5_1_3]=...
            genetic_optimization(u0_generation_5_1_3,u_generation_5_1_3,u_options,...
                round,generation_round,member_number,mutation_probability,exchange_position,parents_number);
        TTS_optimal_5_1_3=TTS_optimal_history_5_1_3(end);
        u_final_5_1_3=(u_final_5_1_3+2)*5;
        disp("TTS of genetic+initiaal points based 30s:");
        disp(TTS_optimal_5_1_3);


%perpendicular search method
        %algorithms settings
        stopping_delta=1;

        %start from random intial points
        u0_5_2_1=u_options(randi(length(u_options))); % k=0
        u_5_2_1=u_options(randi(length(u_options),[1,round]));
        [u0_optimal_5_2_1,u_final_5_2_1,TTS_optimal_history_5_2_1]=perpendicular_search(u0_5_2_1,u_5_2_1,u_options,stopping_delta,round);
        TTS_optimal_5_2_1=TTS_optimal_history_5_2_1(end);
        disp("TTS of perpendicular+initiaal points based on random initial points 2");
        disp(TTS_optimal_5_2_1);

        %start from 30s
        u0_5_2_2=4;
        u_5_2_2=30*ones(1,60);
        [u0_optimal_5_2_2,u_final_5_2_2,TTS_optimal_history_5_2_2]=perpendicular_search(u0_5_2_2,u_5_2_2,u_options,stopping_delta,round);
        TTS_optimal_5_2_2=TTS_optimal_history_5_2_2(end);
        disp("TTS of perpendicular+initiaal points based on initial points 30s");
        disp(TTS_optimal_5_2_2);


    [TTS_5_1,y_5_1,x_5_1,z_5_1,m_5_1,s_5_1,u_5_1] = TTS_calculate([u_final_5_1_3(1),u_final_5_1_3]);
    [TTS_5_2,y_5_2,x_5_2,z_5_2,m_5_2,s_5_2,u_5_2] = TTS_calculate([u_final_5_2_2(1),u_final_5_2_2]); 

    figure('Name','Task 5 traffic situation genetic')
        subplot(2,2,1)
        bar([0:1:round],x_5_1(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_5_1(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_5_1(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_5_1(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")
    
    figure('Name','Task 5 traffic situation perpendicular')
        subplot(2,2,1)
        bar([0:1:round],x_5_2(1,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(u,d)")

        subplot(2,2,3)
        bar([0:1:round],x_5_2(2,:),'b');
        hold on;
        grid on;
        title("number of vehicles,link(o1,d)")

        subplot(2,2,2)
        bar([0:1:round],s_5_2(3,:),'r');
        hold on;
        grid on;
        title("queue length,link(u,d)")

        subplot(2,2,4)
        bar([0:1:round],s_5_2(4,:),'r');
        hold on;
        grid on;
        title("queue length,link(o1,d)")

    figure('Name','Task 5 Genetic Problem')
        subplot(2,3,1)
        plot(TTS_optimal_history_5_1_1)
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_optimal_history_5_1_1)),'r');
        title("random inital point")
        legend('optimized','u(t)=30s')

        subplot(2,3,4)
        bar([0:1:round],[4,u_final_5_1_1])
        title("random inital point:u(t)")
        

        subplot(2,3,2)
        plot(TTS_optimal_history_5_1_2)
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_optimal_history_5_1_2)),'r');
        title("random inital point 2")
        legend('optimized','u(t)=30s')

        subplot(2,3,5)
        bar([0:1:round],[4,u_final_5_1_2])
        title("random inital point 2: u(t)")

        subplot(2,3,3)
        plot(TTS_optimal_history_5_1_3)
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_optimal_history_5_1_3)),'r');
        title("initial points based on 30")
        legend('optimized','u(t)=30s')

        subplot(2,3,6)
        bar([0:1:round],[4,u_final_5_1_3])
        title("initial points based on 30: u(t)")

    figure('Name','Task 5 perpendicular search method')

        subplot(2,2,1)
        plot(TTS_optimal_history_5_2_1)
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_optimal_history_5_2_1)),'r');
        title("random inital point")
        legend('optimized','u(t)=30s')

        subplot(2,2,3)
        bar([0:1:round],[u0_optimal_5_2_1, u_final_5_2_1]);
        title("random inital point: u(t)")

        subplot(2,2,2)
        plot(TTS_optimal_history_5_2_2)
        grid on;
        hold on;
        plot(TTS_ref*ones(1,length(TTS_optimal_history_5_2_2)),'r');
        title(" 30s inital point")
        legend('optimized','u(t)=30s')

        subplot(2,2,4)
        bar([0:1:round],[u0_optimal_5_2_2, u_final_5_2_2]);
        title("30s inital point: u(t)")
        

%% function TTS_calculate in order to calculate TTS

function [TTS,y,x,z,m,s,u] = TTS_calculate(input)
    % calculate the sum of TTS based on input u
    u0=input(1);
    u=input(2:end);
    % u_options=[15,20,25,30,35,40,45];
    % u0=u_options(u0);
    % u=u_options(u);
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
    u=input;
    
    if nargout==1
        TTS;
    end
    

end

function [TTS,y,x,z,m,s,u] = TTS_calculate_2(u0,u)
% calculate the sum of TTS based on input u
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

end

% functions for golden section method
    %golden section is to calculate the local minimum within a given interval, the output can be modified input can corresponding TTS.
function [unew, TTS_current] = golden_section_once(k,u0)

    step = 2*rand()-1;    
    while ((u0(k) + step) >= 45) && ((u0(k) + step) <= 15)
        step = 2*rand()-1;        
    end
  
    lb = u0(k);
    ub = u0(k) + step;
    
    tol = 1.0e-6;
    % if nargin < 4; tol = 1.0e-6; end
    
    nIter = ceil(-2.078087*log(tol/abs(lb-ub)));
  
    golden1 = 0.618033989;
    golden2 = 1 - golden1;
    % First telescoping
    x1 = golden1*ub + golden2*lb;
    x2 = golden2*ub + golden1*lb;
    
    u_1 = u0; u_1(1,k) = x1;
    u_2 = u0; u_2(1,k) = x2;
    u_lb = u0;u_lb(1,k) = lb;
    u_ub = u0;u_ub(1,k) = ub;
    
    f1 = TTS_calculate_2(u_1(1,1), u_1(1,2:61));
    f2 = TTS_calculate_2(u_2(1,1), u_1(1,2:61));
    flb = TTS_calculate_2(u_lb(1,1), u_lb(1,2:61));
    fub = TTS_calculate_2(u_ub(1,1), u_ub(1,2:61));
    
    % 
    for i =1:1:nIter
        if f1 > f2
            lb = x1; x1 = x2; f1 = f2;
            x2 = golden2*lb + golden1*ub;
            u_2(1,k) = x2;
            f2 = TTS_calculate(u_2(1,1), u_1(1,2:61));
            unew = u_2;
            TTS_current = f2;
        elseif f1<f2
            ub = x2; x2 = x1; f2 = f1;
            x1 = golden1*lb + golden2*ub;
            u_1(1,k) = x1;
            f1 = TTS_calculate(u_2(1,1), u_1(1,2:61));
            unew = u_1;
            TTS_current = f1;
        else
            unew = u_ub;
            TTS_current = fub;
        end
    end    
end
%we must be carfeul when dealing with this function's stop certerion, we should count how many times we observe the same output
%when the number is large enough, we stop the searching. 
%haven't finished yet
function [u_output,TTS_current,TTS_optimal_history] = goldensection(u_input,TTS_previous0,golden_iter)
    count_number = 0;
    TTS_previous = TTS_previous0;
    TTS_previous_round = TTS_previous;
    TTS_optimal_history=[];
    while count_number < golden_iter
        for i = 1:1:61
            [u_output, TTS_current] = golden_section_once(i,u_input);
            if(TTS_current > TTS_previous)
            else
                u_input = u_output;
                TTS_previous = TTS_current;
            end 
            TTS_optimal_history=[TTS_optimal_history,TTS_current];
        end
        if TTS_current == TTS_previous_round
            count_number = count_number + 1;
        end
        TTS_previous_round = TTS_current;
    end
end
    
%% function for genetic optimization method
function [u_optimal,TTS_min_history]=genetic_optimization(u0_generation,u_generation,u_options,...
    round,generation_round,member_number,mutation_probability,exchange_position,parents_number)
    %u_optimal for the optimal u(k),TTS_min_history for the TTS history

    u=ones(1,round);
    TTS_generation=ones(member_number,1);
    TTS_min=[];
    for i=1:1:generation_round
        for j=1:1:member_number
            u0=u_options(u0_generation(j));
            for k=1:1:round
                u(k)=u_options(u_generation(j,k));
            end
            TTS(j) = TTS_calculate([u0,u]);
        end
        %sort TTS to find better parents
        sorted_TTS=sort(TTS);
        
        %specified better parents
        parents_position=find(TTS<=sorted_TTS(parents_number),parents_number);
        
        %calculate weight of parents
        w=zeros(1,member_number);
        w(parents_position)=sum(TTS(parents_position))./TTS(parents_position);
        w=w/sum(w);
        
        TTS_min=[TTS_min,min(TTS)];
        pos=find(TTS==min(TTS));
        u_optimal=u_generation(pos(1),:);
        
        %generation children
        children=[];
        for t=1:1:member_number/2
            parent_1=randsample([1:1:member_number],1,true,w);
            parent_2=randsample([1:1:member_number],1,true,w);
            
            %exchange
            child_1=[u_generation(parent_1,1:exchange_position),...
                u_generation(parent_2,exchange_position+1:end)];
            child_2=[u_generation(parent_2,1:exchange_position),...
                u_generation(parent_1,exchange_position+1:end)];
            
            %mutation
            for b=1:1:length(round)
                prob_1=rand(1);
                prob_2=rand(2);
                if (prob_1<=mutation_probability)
                    child_1(b)=randi(length(u_options));
                end
                if (prob_2<=mutation_probability)
                    child_2(b)=randi(length(u_options));
                end
            end
            
            %%add children to new generation
            children=[children;child_1;child_2];
                 
        end
        u_generation=children;
    end
    u_generation;
    TTS_min_history=TTS_min;
end
%% function for perpendicular search method
function [u0_optimal,u_optimal,TTS_min_history]=perpendicular_search(u0,u,u_options,stopping_delta,round)
    %iteration
    delta=2;
%     temp_u0=0;
%     temp_u=zeros(1,round)
    u0_optimal=u0;
    u_optimal=u;
    TTS=[];

    TTS_optimal = TTS_calculate([u0,u]);
    while (delta>=stopping_delta)
    % for r=1:1:5
        for k = 1:1:round+1
            %choose perpendicular search direction
            temp_u=u_optimal;
            direction(k)=1;
            
            %choose the step
            if k==1
               for i= 1:1:length(u_options)
                    temp_u0=u_options(i);
                    [temp_TTS] = TTS_calculate([temp_u0,temp_u]);
                    
                    if temp_TTS <TTS_optimal
                        u0_optimal=temp_u0;
                        TTS_optimal=temp_TTS;
                    end        
                end         
                TTS=[TTS,TTS_optimal];
            else
                for i= 1:1:length(u_options)
                    temp_u(k-1)=u_options(i);
                    temp_TTS = TTS_calculate([temp_u0,temp_u]);
                    if temp_TTS <TTS_optimal
                        u_optimal=temp_u;
                        TTS_optimal=temp_TTS;
                    end    
                end
                TTS=[TTS,TTS_optimal];
            end
            direction(k)=0;
        end
        delta=abs(TTS(end)-TTS(end-round));
    end
    TTS_min_history=TTS;
end

