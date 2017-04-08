%{ 
This is a spring-mass-damper system simulation
Differential equation : my'' + cy' + ky' = u
Transfer function : G(s)= Y(s)/U(s)= 1 / (ms^2 + cs + k)
State-space :
x  = [x1 ; x2]
x' = Ax + Bu
y  = Cx + Du

This code will show
1.spring-mass-damper animation graph
2.x(t),x_dot(t) plot
3.energy plot (kinetic, potential, damping losing, etc.)
%} 


%Restart
clear; % removes all variables from the workspace.
close; % deletes the current figure or the specified figure(s)
%clc;   % clears all input and output from the Command Window display


%Variable/Condition/State-space
m=1;
k=5;
c=4;
Spring_R=0.5;
t=linspace(0,10,100);  % you can change n=100 to n=500 or more to improve UcI accuracy
x0=4;
x0_dot=2;
u=2+0*t;
A=[0 1;-k/m -c/m];
B=[0;1/m];
C=[1 0;0 1];
D=[];
sys=ss(A,B,C,D);
[y,t]=lsim(sys,u,t,[x0;x0_dot]);


%Plot main figure
figure(1);
title('MassDamperSpring');    
for i=1:length(t)    
    %plot X verus t
    subplot(3,1,2);
    plot(y(:,1),t);
    set(gca,'ydir','reverse');
    axis([-5 5 0 10]);
    p0=[y(i,1) t(i)];
    p0_circ=viscircles(p0,0.01,'color','r');
    xlabel('x(m)');
    ylabel('t(s)');
    
    %plot X verus t
    subplot(3,1,3);
    plot(y(:,2),t);
    set(gca,'ydir','reverse');
    axis([-4 4 0 10]);
    p1=[y(i,2) t(i)];
    p1_circ=viscircles(p1,0.01,'color','r');
    xlabel('x_d_o_t(m/s)');
    ylabel('t(s)');
    
    %plot animation
    subplot(3,1,1);
    axis([-5 5 -2 2]);
    %set(gca,'xtick',[],'ytick',[]); %delete the info of axis
    
    %plot animation- ground
    ground_1=line([-5 -5],[-2 2],'color','k','linewidth',3);
    
    %plot animation- SquareMotion
    p1=[-0.5+y(i,1) -1];
    p2=[0.5+y(i,1) -1];
    p3=[0.5+y(i,1) 1];
    p4=[-0.5+y(i,1) 1];
    bottom=line([p1(1) p2(1)],[p1(2) p2(2)]);
    right=line([p2(1) p3(1)],[p2(2) p3(2)]);
    left=line([p3(1) p4(1)],[p3(2) p4(2)]);
    top=line([p4(1) p1(1)],[p4(2) p1(2)]);    

    %plot animation- SpringMotion
    p5=[-0.5+y(i,1) 0.5];
    p6=[(9*-3+(-2.2+y(i,1)))/10 0.5+sqrt(Spring_R.^2-((9*-3+-2.2+y(i,1))/10+3).^2)];
    p7=[(7*-3+3*(-2.2+y(i,1)))/10 0.5-sqrt(Spring_R.^2-(p6(1)+3).^2)];
    p8=[(5*-3+5*(-2.2+y(i,1)))/10 0.5+sqrt(Spring_R.^2-(p6(1)+3).^2)];
    p9=[(3*-3+7*(-2.2+y(i,1)))/10 0.5-sqrt(Spring_R.^2-(p6(1)+3).^2)];
    p10=[(-3+9*(-2.2+y(i,1)))/10 0.5+sqrt(Spring_R.^2-(p6(1)+3).^2)];
    
    spring_1=line([-5 -3],[0.5 0.5],'color','k'); %from ground
    spring_3=line([-3 p6(1)],[0.5 p6(2)],'color','k');
    spring_4=line([p6(1) p7(1)],[p6(2) p7(2)],'color','k');
    spring_5=line([p7(1) p8(1)],[p7(2) p8(2)],'color','k');
    spring_6=line([p8(1) p9(1)],[p8(2) p9(2)],'color','k');
    spring_7=line([p9(1) p10(1)],[p9(2) p10(2)],'color','k');
    spring_8=line([p10(1) -2.2+y(i,1)],[p10(2) p5(2)],'color','k');
    spring_2=line([-2.2+y(i,1) p5(1)],[0.5 p5(2)],'color','k'); %to block    
    
    %plot animation- DamperMotion
    p11=[-0.5+y(i,1) -0.5];
    damper_1=line([-5 -0.5],[-0.5 -0.5],'color','r');    %fixed
    damper_2=line([-0.5 -0.5],[-0.3 -0.7],'color','r');  %fixed
    damper_3=line([-0.7+y(i,1) p11(1)],[-0.5 p11(2)],'color','r');
    damper_4=line([-0.7+y(i,1) -0.7+y(i,1)],[-0.2 -0.8],'color','r'); 
    damper_5=line([-4.7+y(i,1) -0.7+y(i,1)],[-0.2 -0.2],'color','r'); %floor
    damper_6=line([-4.7+y(i,1) -0.7+y(i,1)],[-0.8 -0.8],'color','r'); %bottom
    pause(0.001);
    if(i~=length(t))
        delete(bottom);
        delete(right);
        delete(left);
        delete(top);
        delete(spring_2);
        delete(spring_3);
        delete(spring_4);
        delete(spring_5);
        delete(spring_6);
        delete(spring_7);
        delete(spring_8);
        delete(damper_3);
        delete(damper_4);
        delete(damper_5);
        delete(damper_6);
    end
end
hold off;

%Energy function & calculation
Tk=0.5*k*(y(:,1)).^2;   %Spring Potential Energy
V=0.5*m*(y(:,2)).^2;    %Kinetic energy
UcI=zeros(1,length(t)); %Damping causing energy by trapezoidal integration
Uc=zeros(1,length(t));  %Damping causing energy
Et=zeros(1,length(t));  

for i=1:length(t)
    % calculate damping causing energy by trapezoidal integration
     if(i==length(t))
         UcI(i)=UcI(i-1);
     elseif(i==1)
         UcI(i)=0+c*(t(i+1)-t(i))*((y(i,2)+y(i+1,2))/2).^2;
     else
         UcI(i)=UcI(i-1)+c*(t(i+1)-t(i))*(y(i,2).^2+y(i+1,2).^2)/2;
     end
    % calculate damping causing energy by subtraction
    Et(i)=Tk(i,1)+V(i,1);
    Uc(i)=Et(1)+u(i)*(y(i,1)-4)-Et(i);
end

%Plot energy figure 
figure(2);
xlabel('t(s)');ylabel('energy(J)');title('energy plot');
hold on;
plot(t,Tk);
plot(t,V);
plot(t,Uc(1,:));
plot(t,UcI(1,:));
plot(t,Et(1,:));
legend('Tk','V','Uc','Uc Integration','Etotal');
axis([0 10 0 45]);

for i=1:length(t)
    p0=[t(i) Tk(i,1)];
    p1=[t(i) V(i,1)];
    p2=[t(i) Uc(i)];
    p3=[t(i) Et(i)];
    Tk_circ=viscircles(p0,0.05,'color','k');
    V_circ=viscircles(p1,0.05,'color','k');
    Uc_circ=viscircles(p2,0.05,'color','k');
    Et_circ=viscircles(p3,0.05,'color','k');
    pause(0.001);
    delete(Tk_circ);
    delete(V_circ);
    delete(Uc_circ);
    delete(Et_circ);
end