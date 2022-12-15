%% * Robot (kinematic) model parameters * 
clear all;
close all;

%ta mege8h dinontai se cm
l(1) = 12.0; %to antistoixo l2
l(2) = 11.0; %to antistoixo l4
l(3) = 10.0; %to antistoixo l5

% *sampling period* dinetai se sec
dt = 0.1;

%*DESIRED MOTION PROFILE - TASK SPACE* 
Tf=10.0; %10sec duration of motion 
t=0:dt:Tf;

%start-stop positions in cms
xA = 3.0;	
xB =  9.0; 
yA = 6.00; 
yB = 12.00;
h = 10.00;
zA = h;
zB = h;

%desired trajectory : linear segment (xA,yA,zA)-->(xB,yB,zB)
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); 
disp(' ');   
x(1) = xA; 
y(1) = yA; 

[x,xd,xdd]=jtraj(xA, xB, t);

lambda_x = (xB-xA)/Tf; 
lambda_y = (yB-yA)/Tf; 
a=lambda_y/lambda_x;
disp('Thesis x in certain time');
disp(x);

y = yA + a*(x-xA);
disp('Thesis y in certain time');
disp(y);

yd= a*xd;
disp('Velocity in certain time');
disp(yd);

z = t*(zB - zA)+ zA; 
disp('Thesis z in certain time');
disp(z);

fig0 = figure;
plot(t,xd); 
ylabel('Velocity'); 
xlabel('t (sec)');

fig1 = figure;  
subplot(3,1,1); 
plot(t,x); 
ylabel('x (cm)'); 
xlabel('t (sec)');  

subplot(3,1,2);
plot(t,y); 
ylabel('y (cm)'); 
xlabel('t (sec)');  

subplot(3,1,3);
plot(x,y); 
ylabel('y (cm)'); 
xlabel('x (cm)');  

%% ****** KINEMATIC SIMULATION - Main loop ****** 
disp('Kinematic Simulation ...'); %% 
disp(' '); %%  

%% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE ***** 

%gia thn ar8rwsh q3
den3 = 2*l(2)*l(3);
r3 = x(:).^2 + y(:).^2 + z(:).^2 ;
num3 = r3 + -l(1)^2 -l(2)^2 -l(3)^2;
%q3down = acos(num3/den3); %elbow down
%q3up = acos(-num3/den3); %elbow up
q(:,3) = acos(-num3./den3); %ONE SOLUTION FOR THE SIMULATION

%figq3 = figure;  
%subplot(2,1,1); 
%plot(t,q3up); 
%ylabel('q3up (angle)'); 
%xlabel('t (sec)');

%subplot(2,1,2); 
%plot(t,q3down); 
%ylabel('q3down (angle)'); 
%xlabel('t (sec)');

figq3 = figure;   
plot(t,q(:,3)); 
ylabel('q(:,3) (angle)'); 
xlabel('t (sec)');

%necessary angles from q3 to q2
%ONE SOLUTION FOR THE SIMULATION
s3 = sin(q(:,3));
c3 = cos(q(:,3));
%s3 = sin(q3up);
%s3 = sin(q3down);
%c3 = cos(q3up);
%c3 = cos(q3down);

%gia thn ar8rwsh q2
den2 = ((l(3)^2)*(s3.^2))+((l(2)+l(3)*c3).^2);
r2 = (l(2)+l(3).*c3).*(sqrt((den2-(y(:).^2))));
num2arnhtiko = l(3)*s3.*y(:)-r2;
num2thetiko = l(3)*s3.*y(:)+r2;
%q2down = acos(num2arnhtiko./den2); %elbow down
%q2up = acos(num2thetiko./den2); %elbow down
q(:,2) = acos(num2arnhtiko./den2); %ONE SOLUTION FOR THE SIMULATION

%figq2 = figure;  
%subplot(2,1,1); 
%plot(t,q2up); 
%ylabel('q2up (angle)'); 
%xlabel('t (sec)');

%subplot(2,1,2); 
%plot(t,q2down); 
%ylabel('q2down (angle)'); 
%xlabel('t (sec)');

figq2 = figure;   
plot(t,q(:,2)); 
ylabel('q(:,2) (angle)'); 
xlabel('t (sec)');

%necessary angles from q2 to q1
%for the sum of q2 + q3
%ONE SOLUTION FOR THE SIMULATION
s23 = sin(q(:,2)+q(:,3));
c23 = cos(q(:,2)+q(:,3));
%s23 = sin(q2up+q3up);
%s23 = sin(q2down+q3down);
%c23 = cos(q2up+q3up);
%c23 = cos(q2down+q3down);

%for the q2
%ONE SOLUTION FOR THE SIMULATION
s2 = sin(q(:,2));
c2 = cos(q(:,2));
%s2 = sin(q2up);
%s2 = sin(q2down);
%c2 = cos(q2up);
%c2 = cos(q2down);

%gia thn ar8rwsh q1
den1 = (l(1)^2)+(l(2).*c2+l(3).*c23).^2;
r1 = l(1).*(sqrt(den1-(z(:).^2)));
num1arnhtiko = -(l(2)*c2+l(3)*c23).*z(:)-r1;
num1thetiko = -(l(2)*c2+l(3)*c23).*z(:)+r1;
%q1down = acos(num1arnhtiko./den1); %elbow down
%q1up = acos(num1thetiko./den1); %elbow down
q(:,1) = acos(num1arnhtiko./den1);

%figq1 = figure;  
%subplot(2,1,1); 
%plot(t,q1up); 
%ylabel('q1up (angle)'); 
%xlabel('t (sec)');

%subplot(2,1,2); 
%plot(t,q1down); 
%ylabel('q1down (angle)'); 
%xlabel('t (sec)');

figq1 = figure;   
plot(t,q(:,1)); 
ylabel('q(:,1) (angle)'); 
xlabel('t (sec)');

%% ***** INVERSE DIFFERENTIAL JACOBIAN***** 

%necessary angles from q1 
%ONE SOLUTION FOR THE SIMULATION
s1 = sin(q(:,1));
c1 = cos(q(:,1));
%s1 = sin(q1up);
%s1 = sin(q1down);
%c1 = cos(q1up);
%c1 = cos(q1down);

%boh8hthikes parametroi 
d1 = (l(2).*c2)+(l(3).*c23);
d2 = (l(2).*s2)+(l(3).*s23);

uxe=xd(:);
uye=yd(:);

%h orizousa
dendet = l(3)*l(2).*d1.*s3;
det = 1./dendet;

%for q1
dq(:,1)=det.*((l(3).*c1.*s3).*uxe);

%for q2
a21 = l(3).*c23.*(l(1).*c1+d1.*s1);
a22 = d1*l(3).*s23;
dq(:,2)=det.*((a21.*uxe)+(a22.*uye));

%for q1
a31 = -d1.*(l(1).*c1+d1.*s1);
a32 = -d1.*d2;
dq(:,3)=det.*((a31.*uxe)+(a32.*uye));

figVelocity = figure;  
subplot(3,1,1); 
plot(t,dq(:,1)); 
ylabel('dq(:,1)'); 
xlabel('t (sec)');

subplot(3,1,2); 
plot(t,dq(:,2)); 
ylabel('dq(:,2)'); 
xlabel('t (sec)');

subplot(3,1,3); 
plot(t,dq(:,3)); 
ylabel('dq(:,1)'); 
xlabel('t (sec)');

%% ***** ANIMATION *****
xd1= (s1.*l(2).*c2) + (l(1).*c1);
yd1=l(2).*s2;
xd2=(s1.*c2.*c3.*l(3))-(l(3).*s1.*s2.*s3+l(2).*s1.*c2+l(1).*c1);
yd2=(l(3).*s2.*c3)+(l(3).*c2.*s3)+(l(2).*s2);
kmax = Tf/dt + 1;
figNEW = figure; 
axis([-20 10 -10 10]) %%set xy plot axes (caution: square axes, i.e. dx=dy) 
axis on 
hold on 
xlabel('x (cm)'); 
ylabel('y (cm)');
plot(xd,yd,'rs'); 
dtk=1000; %% plot robot position every dtk samples, to animate its motion 
plot([0],[0],'o');
for tk=1:dtk:kmax,    %%% 	
   pause(0.1);	%% pause motion to view successive robot configurations    
   plot([0,xd1(tk)],[0,yd1(tk)]);					
   plot([xd1(tk)],[yd1(tk)],'o');    
   plot([xd1(tk),xd2(tk)],[yd1(tk),yd2(tk)]);	
   plot([xd2(tk)],[yd2(tk)],'y*');    
   plot([xd(tk)],[yd(tk)],'g+');  
end  