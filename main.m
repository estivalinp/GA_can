%独立展示矩阵和运算矩阵
%重新检查适应度√
%2015.10.26 23:10
%检查过程√
%循环√
%函数
%交叉、变异后检查父代？
%加入罚函数？√（简单淘汰法）
%2015.10.27 20:16
%2015.10.29 11:47
%添加同步显示Task Manager√
%▲改正交叉算子错误！！！√
%2015.10.29 14:09
%2015.10.29 18:39
%修改图形参数显示
%2015.10.29 22:32
clear all
clc
tic
%Parameter definition
global popu_num;
global code_length;
global genal_num;
global genal_count;
global popu;
global crossover_probability;
global mutation_probability;

global per_best;

%Basic parameter
x_min = 0;
x_max = 31;
accuracy = 0.001;%0.001
cost = 0.0654; % Cost of can material per cm^2

popu_num = 4000;%1000
code_length = ceil(log2( (x_max-x_min) / accuracy ));
genal_num = 100;
genal_count = 1;
popu = zeros(6,popu_num); % Initialize population to zeros.
crossover_probability = 0.5;
mutation_probability = 0.01;

mask1 = 2^(2*code_length)-(2^code_length-1); % Substring mask
mask2 = 2^code_length-1;
shift = code_length;  % Bit shifting

%disp('0.Initialize population(Coding)');
count=1; % Counter of number of population.
while count ~= popu_num + 1
    % The second method.(faster for Matlab)
    popu(1, count) = randi([0 2^(code_length * 2)-1],1,1); % random integer
    % Is an enable solution?
    flag = 0; % Marker for is(able) or is not(disable) make a new string.
    d_1 =  bitshift(bitand(popu(1,count), mask1), -shift);  %Binary.
    h_1 = bitand(popu(1,count), mask2);
    d = x_min + (x_max-x_min)/(2^code_length-1)*d_1; %bin2dec
    h = x_min + (x_max-x_min)/(2^code_length-1)*h_1;
    % Do not have 300 ml volume?
    if (pi * d^2 * h)/4 < 300
%         %disp('<300');
        flag = 1;
    end
    % Have 300 ml volume.
    if flag == 0    
        % Repeat? && Zero?
        if count==1 && ( bitand(popu(1,count), mask2)==0 || ... % first one
            bitand(popu(1,count), mask1)==0 );
            flag = 1;
%             %disp('!');
        end
        for k = 1:count-1
            if popu(1,count) == popu(1,k) || ... % next one
                    bitand(popu(1,count), mask2)==0 ||...
                    bitand(popu(1,count), mask1)==0   
                flag = 1;
%                 %disp('!');
                break
            end
        end
    end
    
    if flag == 0 % Despair
        % %disp(dec2bin(popu(2,count),code_length*2))
        count = count + 1;
    else % Rebuid
        popu(1,count) = 0;
    end
    % %disp(' ');
end
%disp('-----------------------------------');

while genal_count <= genal_num
    %disp('1.Decode and caculate fitness value');    
    for i = 1:popu_num    
        popu(3,i) = bitshift(bitand(popu(1,i), mask1), -shift); %d
        popu(4,i) = bitand(popu(1,i), mask2); %h
        popu(3,i) = x_min + (x_max-x_min)/(2^code_length-1)*popu(3,i);% d 
        popu(4,i) = x_min + (x_max-x_min)/(2^code_length-1)*popu(4,i);% h
        popu(6,i) = (pi * popu(3,i)^2 * popu(4,i))/4; % Volume 
        if popu(6,i) < 300
            popu(5,i) = inf; % penalized
        else            
        	popu(5,i) = cost * (pi*popu(3,i)^2/2 + pi*popu(3,i)*popu(4,i)); % fitness value   
        end
        % popu(5,i) = cost * (pi*popu(3,i)^2/2 + pi*popu(3,i)*popu(4,i)); % fitness value    
        % popu(6,i) = (pi * popu(3,i)^2 * popu(4,i))/4; % Volume    
    end        
    %disp('Initial population');
    %disp(dec2bin(popu(1,:), code_length*2));   
    %disp('Initial average fitness')
    %disp(mean(popu(5,:)));
    %disp('-----------------------------------');
    
    %Draw--------------------------------------
    per_best(1,genal_count) = min(popu(5,:));
    per_best(2,genal_count) = mean(popu(5,popu(5,:)~=inf));  
    mutation_point = find(popu(2,:)==1);
    best_point = find(popu(5,:)==per_best(1,genal_count));
    % Figura1 (GA Population Point)
    subplot(1,3,1)
    plot(popu(3,:),popu(4,:),'.','MarkerSize',4)
    hold on
    plot(popu(3,mutation_point),popu(4,mutation_point),'.','MarkerSize',4,'color','red')
    hold on
    plot(popu(3,best_point),popu(4,best_point),'.','MarkerSize',4,'color','green')
    text(20,28,['generation: ',num2str(genal_count)],'FontWeight','bold');
    text(20,26,['best value: ',num2str(per_best(1,genal_count))],'FontWeight','bold');
    plot([19,19,29,29,19],[25,29,29,25,25],'black');
    title('GA Population')
    xlabel('D')
    ylabel('H')
    axis([0 31 0 31]);
    set(gca,'DataAspectRatio',[1 1 1])
    hold off
    if genal_count > 1
        % Figura2 (GA Optimization Curves)
        subplot(1,3,2)  
        hold on
        %grid on
        plot([genal_count-2,genal_count-1], [per_best(1,genal_count-1),per_best(1,genal_count)],'red')   
        hold on        
        plot([genal_count-2,genal_count-1], [per_best(2,genal_count-1),per_best(2,genal_count)]) 
        axis([0 100 10 110]);
        set(gca,'DataAspectRatio',[1 1 1])
        title('GA Optimization Curves')
        xlabel('Generation')
        ylabel('Cost')
        % Figura3 (GA Convergence Curves)
        subplot(1,3,3)
        hold on        
        plot([genal_count-2,genal_count-1], [per_best(1,genal_count-1),per_best(1,genal_count)],'red') 
        axis([0 100 16.2 17]);
        set(gca,'DataAspectRatio',[125 1 1])
        title('GA Convergence Curves')
        xlabel('Generation')
        ylabel('Cost')
    end
    pause(0.1)
    hold off
    %Draw--------------------------------------

    %disp('2.Reproduction Operator');
    popu(:,:) = popu(:,randperm(numel(popu(1,:)))); % Randomly arranged.
    %disp('Mixed');
    %disp(dec2bin(popu(1,:), code_length*2));
    %disp(' ');
    popu_tmp = popu;
    % Compare fitness(copy or eliminate)using tournament selection.
    for i = 1:popu_num
        if popu_tmp(5,i) < popu_tmp(5,mod(i,popu_num)+1) % 2 3 4 5 6 1
            popu(:,i) = popu_tmp(:,i);
        else
            popu(:,i) = popu_tmp(:,mod(i,popu_num)+1);
        end
    end
    %disp('After');
    %disp(dec2bin(popu(1,:), code_length*2));
    %disp('-----------------------------------');

    %disp('3.Crossover Operator');
    popu(:,:) = popu(:,randperm(numel(popu(1,:)))); % Randomly arranged.
    %disp('Mixed');
    %disp(dec2bin(popu(1,:), code_length*2));
    %disp(' ');
    for i=1:floor(popu_num/2)
        if rand() < crossover_probability
            cross_point = randi([1 code_length*2-1]); % Get crossover point ramdonly
            substring1 = bitand(popu(1,i), 2^(2*code_length-cross_point)-1);
            substring2 = bitand(popu(1,popu_num+1-i), 2^(2*code_length-cross_point)-1);
            popu(1,i) = popu(1,i) - substring1 + substring2; % exchange
            popu(1,popu_num+1-i) = popu(1,popu_num+1-i)-substring2+substring1;
        end
    end
    %disp('After');
    %disp(dec2bin(popu(1,:), code_length*2));
    %disp('-----------------------------------');

    %disp('4.Mutation Operator');
    mutation_array = rand(popu_num, code_length*2)<mutation_probability; % random float
    for i=1:popu_num    
        %disp(dec2bin(popu(1,i), code_length*2));
        popu_temp = popu(1,i);
        popu(2,i) = 0; % No mutation fistly.
        popu(1,i) = bitxor(popu(1,i),bi2de(mutation_array(i,:)));%!!!!!!!!!!!!!!!!!!
        if(popu_temp-popu(1,i)~=0)
            popu(2,i) = 1; % Mark the mutation point.
        end
        %disp(dec2bin(popu(1,i), code_length*2));
%          if bi2de(mutation_array(i,:)) ~= 0
%             %disp('^');
%         else
%             %disp(' ');
%          end
    end
    %disp('-----------------------------------');   
    
    genal_count = genal_count+1;
end

format bank
%disp('Best solution');
best_side = min(popu(5,:));
%disp(best_side(1));
popu=popu';

toc

