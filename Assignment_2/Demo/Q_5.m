%constant
c=60;
round=60;
l_veh=7;

E1=9;
E2=13;
E3=5;

%discrete green light time set

u_options=[15,20,25,30,35,40,45];



%algorithm settins

generation_round=50;

member_number=50;

%the exchange_position of parents
exchange_position=30;

%mutation probability
mutation_probability=0.5;

%initial generation 0

u_generation=randi(length(u_options),[1,round]);
u0_generation=randi(length(u_options),[member_number,1])
for i =2:1:round
    u_generation=[u_generation;randi(length(u_options),[1,round])];
end



%iteration
u=ones(1,round);
y_generation=ones(member_number,1);
for i=1:1:2
    for j=1:1:member_number
        u0=u_options(u0_generation(j));
        for k=1:1:round
            u(k)=u_options(u_generation(j,k));
        end
        y(j)=TTS_calculate(u0,u);
    end
    %sort y to find better parents
    sorted_y=sort(y);
    
    %specified better parents
    parents_position=find(y>=sorted_y(20),20);
    
    %generation children
    children=[];
    for t=1:1:(member_number-20)/15
        parent_1=parents_position(randi(length(parents_position)));
        parent_2=parents_position(randi(length(parents_position)));
        
        %exchange
        child_1=[u_generation(parent_1,1:exchange_position),...
            u_generation(parent_2,exchange_position+1:end)];
        child_2=[u_generation(parent_2,1:exchange_position),...
            u_generation(parent_1,exchange_position+1:end)];
        
        %mutation
        for b=1:1:length(round)
            prob_1=rand(1);
            prob_2=rand(2);
            if (prob_1>=mutation_probability)
                child_1(b)=randi(length(u_options));
            end
            if (prob_2>=mutation_probability)
                child_2(b)=randi(length(u_options));
            end
        end
        
        %%add children to new generation
        children=[children;child_1;child_2]
             
    end
    size(children)
    size(u_generation(parents_position,:))
    u_generation=[u_generation(parents_position);children];
    
end
