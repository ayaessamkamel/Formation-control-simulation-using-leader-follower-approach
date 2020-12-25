

dt=0.01;%sampling period

iteration =100;%numbers of iterations 

z0 = [5 18 8 7 10 2]'%initial condition
%desired relative position
% c12 = [-1,-1];
% c13 = [1,-1];
% c23 = [2,0];
% c32 = [-2,0];
c12 = [1,1.5];
c13 = [-1,1.5];
c23 = [-2,0];
c32 = [2,0];

%desired vectors of relative positions (c_mat = [sum of the desired vectors of positions for x1, y1, x2, y2, x3, y3])
c_mat = [0,0,c12(1,1) + c32(1,1),c12(1,2) + c32(1,2),c13(1,1) + c23(1,1), c13(1,2) + c23(1,2)]';

%initialize the array for the state of agents
%z = [z1x(0) z1x(k) ...
%       z1y(0) z1y(k) ......
%size(z) = [numNodes*2, k]
z = zeros(numNodes*2, iteration+1);

%Initialize the array for the distance between each 'adjacent' agents
% dis = [  0   dis12 dis13 ;
%        dis21   0   dis23 ;
%        dis31 dis32   0    ...];
dis = zeros(numNodes,numNodes,length(z));

%Initialize the array for the RPF's derivative
% d_v = [  0   dv_12 dv_13 ;
%        dv_21   0   dv_23 ;
%        dv_31 dv_32   0    ...];
d_v =zeros(numNodes,numNodes,length(z));
z(:,1)=z0; 
for k =1:iteration+1
    
    dis(1,2,k) = norm(z(1:2,k) - z(3:4,k)); %distance between node 1 and node 2
    dis(1,3,k) = norm(z(1:2,k) - z(5:6,k)); %distance between node 1 and node 3
    dis(2,1,k) = dis(1,2,k);                    %distance between node 2 and node 1 (same as dis21)
    dis(2,3,k) = norm(z(3:4,k) - z(5:6,k)); %distance between node 2 and node 3
    dis(3,1,k) = dis(1,3,k);                    %distance between node 3 and node 1
    dis(3,2,k) = dis(2,3,k);                    %distance between node 3 and node 2
    
    d_v(1,2,k) = -4*eta*(d^2 - dis(1,2,k)^2) / (d^2*dis(1,2,k)^5); %potential derivative of node 1 to 2
    d_v(2,1,k) = -d_v(2,1,k);                                      %potential derivative of node 2 to 1
    d_v(1,3,k) = -4*eta*(d^2 - dis(1,3,k)^2) / (d^2*dis(1,3,k)^5); %potential derivative of node 1 to 3
    d_v(3,1,k) = -d_v(1,3,k);
   
    d_v(2,3,k) = -4*eta*(d^2 - dis(2,3,k)^2) / (d^2*dis(2,3,k)^5);
    d_v(3,2,k) = -d_v(2,3,k);
   
    % if-else condition, make sure that RPF starts to activate when obstacles enter the danger zone
    % if (current distance)^2 - d^2 <=0
    %    de = 1 
    % else
    %    de = 0
    de12 = (dis(1,2,k)^2 - d^2 <= 0);
    de13 = (dis(1,3,k)^2 - d^2 <= 0);
    de23 = (dis(2,3,k)^2 - d^2 <= 0);
    
    
    % second stage: move together in certain pattern
    if ((z(3,k) - z(1,k)) + (z(3,k) - z(5,k)) - c_mat(3,1) < 0.5 ...
            && (z(4,k) - z(2,k) + (z(4,k) - z(6,k)) - c_mat(4,1) <0.5)...
            && (z(5,k) - z(1,k) + (z(5,k) - z(3,k)) - c_mat(5,1) < 0.5)...
            && (z(6,k) - z(2,k) + (z(6,k) - z(4,k)) - c_mat(6,1) < 0.5))
        % leader bot moving condition
        if (k <= 100)
            k_x = 1;
            k_y = 1.005;
        elseif (k > 100 && k <= 150)
            k_x = 1.005;
            k_y = 1.005;
        elseif (k > 150 && k <= 200)
            k_x = 1;
            k_y = 1.005;
        elseif (k > 200 && k <= 250)
            k_x = 0.995;
            k_y = 1.005;
        elseif (k > 250 && k <= 300)
            k_x = 1.005;
            k_y = 1.005;
        elseif (k > 300)
            k_x = 1;
            k_y = 1.005;           
        end
    else
        k_x = 1;
        k_y = 1.001;  
    end
  
    %nonlinear equation of the agent dynamic (z(k+1) = z(k) - APF's*dt - RPF's*dt))
    z(1,k+1) = k_x * z(1,k); 
    
    z(2,k+1) = k_y * z(2,k);
    
    z(3,k+1) = z(3,k) - de12*d_v(2,1,k)*dt - de23*d_v(2,3,k)*dt...
               - kd*( (z(3,k) - z(1,k)) + (z(3,k) - z(5,k)) - c_mat(3,1) )*dt; %z2_x
    
    z(4,k+1) = z(4,k) - de12*d_v(2,1,k)*dt - de23*d_v(2,3,k)*dt...
               - kd*( (z(4,k) - z(2,k)) + (z(4,k) - z(6,k)) - c_mat(4,1) )*dt; %z2_y
    
    z(5,k+1) = z(5,k) - de13*d_v(3,1,k)*dt - de23*d_v(3,2,k)*dt...
               - kd*( (z(5,k) - z(1,k)) + (z(5,k) - z(3,k)) - c_mat(5,1) )*dt; %z3_x
    
    z(6,k+1) = z(6,k) - de13*d_v(3,1,k)*dt - de23*d_v(3,2,k)*dt...
               - kd*( (z(6,k) - z(2,k)) + (z(6,k) - z(4,k)) - c_mat(6,1) )*dt; %z3_y
    

    
end

figure(2)
plot(z(1,1),z(2,1),'ro',z(3,1),z(4,1),'bo',z(5,1),z(6,1),'go')
hold on
for i=2:iteration
    pause(0.005);
    fprintf('%d\n', i);
    plot(z(1,i),z(2,i),'r:.')
    hold on
    plot(z(3,i),z(4,i),'b:.')
    hold on
    plot(z(5,i),z(6,i),'g:.')
    drawnow
    
    hold on
    grid on
    xlabel('x')
    ylabel('y')
    grid on
end
plot(z(1,iteration+1), z(2,iteration+1),'xr',z(3,iteration+1), z(4,iteration+1),'xb',...
     z(5,iteration+1), z(6,iteration+1),'xg');
legend('Agent1','Agent2','Agent3')
