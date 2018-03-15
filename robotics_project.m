%% Directs Kinematics
clear all;
close all;
clc;

syms theta1 theta2 theta3 d4;

c1 = cos(theta1);
c2 = cos(theta2);
c3 = cos(theta3);

s1 = sin(theta1);
s2 = sin(theta2);
s3 = sin(theta3);
% Matrix Malculation
A01 = [c1 0 -s1 0; s1 0 c1 0; 0 -1 0 0; 0 0 0 1];
A12 = [c2 0 s2 0; s2 0 -c2 0; 0 1 0 0; 0 0 0 1];
A23 = [c3 -s3 0 0; s3 c3 0 0; 0 0 1 0; 0 0 0 1];
A34 = [1 0 0 0; 0 1 0 0; 0 0 1 d4; 0 0 0 1];
 
A02 = A01*A12;
A03 = A02*A23;
A04 = A03*A34;
A13 = A12*A23;
A14 = A13*A34;
A24 = A23*A34;

%% Circle Points

goal = 0:pi/20:2*pi;
theta = pi/6;

x = 8.*cos(goal).*cos(theta);
y = 8.*sin(goal);
z = 8.*cos(goal).*sin(theta);
%% Inverse Kinematics
%υπολογίζουμε την αντίστροφη κινηματική
q = [];
for i=1:length(goal)
    q = [q; atan2(y(i), x(i)) pi/2-atan2(z(i),(x(i)^2+y(i)^2)^0.5) 0 (x(i)^2+y(i)^2+z(i)^2)^0.5];
    if q(i, 1) < 0
        q(i, 1) = 2*pi+q(i, 1);
    end
end
q = q';

coord = [];
%Ο πίνακας coord είναι απλά για επαλήθευση της ορθής κινηματικής, έχει 10 στήλες
%και σε κάθε γραμμή μια συντεταγμένη για τα 10 σημεία της τροχιάς μας



x = 8.*cos(goal).*cos(theta);
y = 8.*sin(goal);
z = 8.*cos(goal).*sin(theta);

for i=1:length(goal)
    
    theta1 = q(1, i);
    theta2 = q(2, i);
    theta3 = q(3, i);
    
    d4 = q(4, i);
    A04k = eval(A04);
    
    coord = [coord A04k*[0 0 0 1]'];   
end

plot3(x,y,z, ':ro');
hold on;
plot3(coord(1, :), coord(2, :), coord(3, :));
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Επιθυμητή τροχιά της αρπάγης');


%% Velocity

close all;

r = [];
%ο r είναι ο πίνακας με τα διανύσματα θεσης κάθε σημείου και ουσιαστικά τον βρίσκουμε 
% για να κάνουμε το εξωτερικό γινόμενο με τη γωνιακή ταχύτητα και να βρίσκουμε κάθε φορά τη γραμμική ταχύτητα

for i=1:length(goal)
    rx = x(i)/8;
    ry = y(i)/8;
    rz = z(i)/8;
    
    r = [r ; rx ry rz];
end

a1 = [x(2)-x(1) y(2)-y(1) z(2)-z(1)];
a1 = a1./norm(a1);

wmega = cross(r(1, :), a1);
wmega = wmega.*2*pi/norm(wmega);

%O πίνακας v περιέχει τις γραμμικές ταχύτητες
%v = [0 0 0];
v= [];

for i=1:length(goal)
    v = [v ; cross(wmega, r(i, :))];
end

%% Jacobian
%υπολογισμός της ιακωβιανής κάθε βαθμού ελευθερίας με εξωτερικό γινόμενο
J1 = [cross([0 0 1]', A04(1:3, 4)); [0 0 1]'];

J2 = [cross(A01(1:3, 3), A04(1:3, 4)-A01(1:3, 4)); A01(1:3, 3)]; 

J3 = [cross(A02(1:3, 3), A04(1:3, 4)-A02(1:3, 4)); A02(1:3, 3)];

J4 = [A03(1:3, 3); [0 0 0]'];
%δημιουργούμε τον τελικό ιακωβιανό πίνακα 
J = [J1 J2 J3 J4];

%% Υπολογισμός ταχύτητας
qdot = zeros(4, length(goal));

for i=1:length(goal)
    
    theta1 = q(1, i);
    theta2 = q(2, i);
    theta3 = q(3, i);
    d4 = q(4, i);
    
    Ja = eval(J);
    Ja = pinv(Ja);
    qdot(:, i) = Ja*[wmega v(i, :)]';   
end

%% Πολυώνυμα παρεμβολής
clc;
zer = zeros(16, length(goal));

time = 0:8/(length(goal)-1):8;

for i=1:4
    zer_row = zeros(4, length(goal));
    for j=1:length(goal)
        if j==length(goal)           
            pol = [time(j)^3 time(j)^2 time(j) 1; 3*time(j)^2 2*time(j) 1 0; time(1)^3 time(1)^2 time(1) 1; 3*time(1)^2 2*time(1) 1 0];
            zer_row(:, j) = pol\[q(i, j) qdot(i, j) q(i, 1) qdot(i, 1)]';
        else
            pol = [time(j)^3 time(j)^2 time(j) 1; 3*time(j)^2 2*time(j) 1 0; time(j+1)^3 time(j+1)^2 time(j+1) 1; 3*time(j+1)^2 2*time(j+1) 1 0];
            zer_row(:, j) = pol\[q(i, j) qdot(i, j) q(i, j+1) qdot(i, j+1)]';
        end
    end
    zer((i-1)*4+1:i*4, :) = zer_row;
end

syms t;

theta_1 = sym(zeros(1, length(goal)));
theta_2 = sym(zeros(1, length(goal)));
theta_3 = sym(zeros(1, length(goal)));
d_4 = sym(zeros(1, length(goal)));

for i = 1:length(goal)
    
    theta_1(i) = zer(1, i)*t^3+zer(2, i)*t^2+zer(3, i)*t+zer(4, i);
    
    theta_2(i) = zer(5, i)*t^3+zer(6, i)*t^2+zer(7, i)*t+zer(8, i);
    
    theta_3(i) = zer(9, i)*t^3+zer(10, i)*t^2+zer(11, i)*t+zer(12, i);
    
    d_4(i) = zer(13, i)*t^3+zer(14, i)*t^2+zer(15, i)*t+zer(16, i);
end


%% Figures
clc;
close all;
P = [];

goal = 0:pi/20:2*pi;
theta = pi/6;

x = 8.*cos(goal).*cos(theta);
y = 8.*sin(goal);
z = 8.*cos(goal).*sin(theta);


for i=1:length(goal)
    theta1 = theta_1(i);
    theta2 = theta_2(i);
    theta3 = theta_3(i);
    d4 = d_4(i);
    P = [P eval(A04)*[0 0 0 1]'];
end

time = 0:8/(length(goal)-1):8;
Pxsimul = [];
Pysimul = [];
Pzsimul = [];

for i =1:length(goal)
    
    t = time(i);
    Pxsimul = [Pxsimul eval(P(1, i))];
    Pysimul = [Pysimul eval(P(2, i))];
    Pzsimul = [Pzsimul eval(P(3, i))];
end


%% Simulation
timer = 0:8/161:8;

Px = [];
Py = [];
Pz = [];
for i =1:length(timer)
    t = timer(i);
    if timer(i) == 0
        Px = [Px eval(P(1, 1))];
        Py = [Py eval(P(2, 1))];
        Pz = [Pz eval(P(3, 1))];
    elseif timer(i) == 8
        Px = [Px eval(P(1, length(goal)))];
        Py = [Py eval(P(2, length(goal)))];
        Pz = [Pz eval(P(3, length(goal)))];
    else
        for j=1:length(goal)       
            if timer(i)>(j-1)*8/(length(goal)-1) && timer(i)<j*8/(length(goal)-1)
                Px = [Px eval(P(1, j))];
                Py = [Py eval(P(2, j))];
                Pz = [Pz eval(P(3, j))];
            end
        end
    end
end

theta__1=[];
theta__2=[];
theta__3=[];
d__4=[];
for i =1:length(timer)
    t = timer(i);
    
    if timer(i) == 0
        theta__1=[theta__1 eval(theta_1(1))];
        theta__2=[theta__2 eval(theta_2(1))];
        theta__3=[theta__3 eval(theta_3(1))];
        d__4=[d__4 eval(d_4(1))];
    elseif timer(i) == 8
        theta__1=[theta__1 eval(theta_1(length(goal)))];
        theta__2=[theta__2 eval(theta_2(length(goal)))];
        theta__3=[theta__3 eval(theta_3(length(goal)))];
        d__4=[d__4 eval(d_4(length(goal)))];
    else
        for j=1:length(goal)        
            if timer(i)>=(j-1)*8/(length(goal)-1) && timer(i)<=j*8/(length(goal)-1)
               theta__1=[theta__1 eval(theta_1(j))];
               theta__2=[theta__2 eval(theta_2(j))];
               theta__3=[theta__3 eval(theta_3(j))];
               d__4=[d__4 eval(d_4(j))];
            end
        end
    end
end


plot(timer, Px, '*-c');
title('γράφημα Px σε σχέση με το χρόνο t');
xlabel('Χρόνος(sec)');
ylabel('Px');
hold on;
plot(time, Pxsimul, 'or');
plot(time, x, '-.k');
figure;

plot(timer, Py, '*-c');
title('γράφημα Py σε σχέση με το χρόνο t');
xlabel('Χρόνος(sec)');
ylabel('Py');
hold on;
plot(time, Pysimul, 'or');
plot(time, y, '-.k');
figure;

plot(timer, Pz, '*-c');
title('γράφημα Pz σε σχέση με το χρόνο t');
xlabel('Χρόνος(sec)');
ylabel('Pz');
hold on;
plot(time, Pzsimul, 'or');
plot(time, z, '-.k');
figure;

plot3(Px, Py, Pz, '*-m');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Επιθυμητή τροχιά αρπάγης');
hold on;
plot3(Pxsimul, Pysimul, Pzsimul, 'ok');
plot3(x, y, z, '-.');
legend('Πραγματική τροχιά','σημεία τροχιάς','επιθυμητή τροχιά','Location','northwest');
figure;

plot(timer, theta__1, 'd-m');
title('θ1(t)');
xlabel('Χρόνος(sec)');
ylabel('θ1');
hold on;
plot(time, q(1,:), '-.k');
legend('θ1(t)', 'σημεία επιθυμητής τροχιάς','Location','southeast');
figure;

plot(timer, theta__2, 'd-m');
title('θ2(t)');
xlabel('Χρόνος(sec)');
ylabel('θ2');
hold on;
plot(time, q(2,:), '-.k');
legend('θ2(t)', 'σημεία επιθυμητής τροχιάς');
figure;

plot(timer, theta__3, 'd-m');
title('θ3(t)');
xlabel('Χρόνος(sec)');
ylabel('θ3');
hold on;
plot(time, q(3,:), '-.k');
legend('θ3(t)', 'σημεία επιθυμητής τροχιάς','Location','southoutside');
figure;

plot(timer, d__4, 'd-m');
title('d4(t)');
xlabel('Χρόνος(sec)');
ylabel('d4');
hold on
plot(time, q(4,:), '-.k');
legend('d(t)', 'σημεία επιθυμητής τροχιάς','Location','southeast');
figure;


%% Dynamics

Qr = [0 -1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0 ]; %Q περιστροφικός
Qp = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0];   %Q πρισματικός

diag = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
A = [A01 A02 A03 A04; diag A12 A13 A14; diag diag A23 A24; diag diag diag A34];
Q = [Qr; Qr; Qr; Qp];
U_ij = [];
U_ijk = [];
% Η πρώτη for (γραμμές 294-307) είναι για τον υπολογισμό των Uij
% (βιβλιογραφία Robotics Control, Sensing, Vision and Intelligence Fu, Gonzalez and Lee

for i=1:4
    test = [];
    for j=1:4
        if j>i
          Uij = zeros(4,4);
        elseif j==1
            Uij = Q((j-1)*4+1:(j-1)*4+4, :)*A(1:4, (i-1)*4+1:(i-1)*4+4);
        else 
            Uij = A(1:4,(j-2)*4+1:(j-2)*4+4)*Q((j-1)*4+1:(j-1)*4+4, :)*A((j-1)*4+1:(j-1)*4+4, (i-1)*4+1:(i-1)*4+4);
        end
        test = [test Uij];
    end
    U_ij = [U_ij; test];
end
% Η δεύτερη for ( γραμμές 310-334) είναι για τον υπολογισμό των Uijk σύμφωνα με τη θεωρία 

for i=1:4
    test = [];
    for j=1:4
        for k=1:4
            if j>i || k>i
              Uijk = zeros(4,4);
            elseif j == 1 && k == 1
                Uijk = Q((j-1)*4+1:(j-1)*4+4, :)*Q((k-1)*4+1:(k-1)*4+4, :)*A(1:4, (i-1)*4+1:(i-1)*4+4);
            elseif j == 1
                Uijk = Q((j-1)*4+1:(j-1)*4+4, :)*A(1:4, (k-2)*4+1:(k-2)*4+4)*Q((k-1)*4+1:(k-1)*4+4, :)*A((k-1)*4+1:(k-1)*4+4, (i-1)*4+1:(i-1)*4+4);
            elseif k == 1
                Uijk = Q((k-1)*4+1:(k-1)*4+4, :)*A(1:4, (j-2)*4+1:(j-2)*4+4)*Q((j-1)*4+1:(j-1)*4+4, :)*A((j-1)*4+1:(j-1)*4+4, (i-1)*4+1:(i-1)*4+4);
            elseif j == k
                Uijk = A(1:4, (j-2)*4+1:(j-2)*4+4)*Q((j-1)*4+1:(j-1)*4+4, :)*Q((k-1)*4+1:(k-1)*4+4, :)*A((k-1)*4+1:(k-1)*4+4, (i-1)*4+1:(i-1)*4+4);
           elseif j > k
                Uijk = A(1:4, (k-2)*4+1:(k-1)*4)*Q((k-1)*4+1:k*4, :)*A((k-1)*4+1:k*4, (j-2)*4+1:(j-1)*4)*Q((j-1)*4+1:j*4, :)*A((j-1)*4+1:j*4, (i-1)*4+1:i*4);
           else
                Uijk = A(1:4, (j-2)*4+1:(j-1)*4)*Q((j-1)*4+1:j*4, :)*A((j-1)*4+1:j*4, (k-2)*4+1:(k-1)*4)*Q((k-1)*4+1:k*4, :)*A((k-1)*4+1:k*4, (i-1)*4+1:i*4);

           end
            test = [test Uijk];
        end
    end
    U_ijk = [U_ijk; test];
end

%% D
I = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1];

D = [];

for i=1:4
    Da = [];
    for k=1:4
        Dik = 0;
        for j=max([i k]):4
            Dik = Dik + trace(U_ij((j-1)*4+1:(j-1)*4+4,(k-1)*4+1:(k-1)*4+4)*I*U_ij((j-1)*4+1:(j-1)*4+4,(i-1)*4+1:(i-1)*4+4)');
        end
        Da = [Da Dik];
    end
    D = [D; Da];
end

%% hikm
h_ikm = [];

for i=1:4
    ha = [];
    for k=1:4
        for m=1:4
            hikm = 0;
            for j=max([i k m]):4
                hikm = hikm + trace(U_ijk((j-1)*4+1:(j-1)*4+4,(k-1)*16+(m-1)*4+1:(k-1)*16+(m-1)*4+4)*I*U_ij((j-1)*4+1:(j-1)*4+4,(i-1)*4+1:(i-1)*4+4)');
            end
            ha = [ha hikm];
        end
    end
    h_ikm = [h_ikm; ha];
end

%% h

qsym = vpa([theta_1; theta_2; theta_3; d_4], 5);

qsymd = vpa(diff(qsym), 5);

qsymdd = vpa(diff(qsymd), 5);

h = [];
for i=1:4
    hi = [];
    for j=1:length(goal)
        hj = 0;
        for k=1:4
            for m=1:4
                hj = hj + h_ikm(i,(k-1)*4+m)*qsymd(k, j)*qsymd(m, j);
            end
        end
        hi = [hi hj];
    end
    h = [h; hi];
end
            
        
%% c
c = [];

g = [0 0 -9.81 0]; %διάνυσμα βαρύτητας

r = [0 0 0 1]'; 

for i=1:4
    ci = 0;
    for j=1:4
        ci = ci - g*U_ij((j-1)*4+1:(j-1)*4+4,(i-1)*4+1:(i-1)*4+4)*r;
    end
    c = [c; ci];
end

%% Torques

torque = [];
torque_th1 = [];
torque_th2 = [];
torque_th3 = [];
torque_d4 = [];

for j=1:length(goal)
	torque = [torque D*qsymdd(:, j)+h(:, j)+c];
end
torque__ = eval(torque);
for i =1:length(goal);
    
    t = time(i);
    torque_th1 = [torque_th1 eval(torque__(1,i))];
    torque_th2 = [torque_th2 eval(torque__(2,i))];
    torque_th3 = [torque_th3 eval(torque__(3,i))];
    torque_d4  = [torque_d4 eval(torque__(4,i))];
end
plot(time, torque_th1(1, :), '*-g');
title('ροπή της θ1');
xlabel('Χρόνος(sec)');
ylabel('Ροπή θ1');
figure;

plot(time, torque_th2(1, :), '*-g');
title('ροπή της θ2');
xlabel('Χρόνος(sec)');
ylabel('Ροπή θ2');
figure;

plot(time, torque_th3(1, :), '*-g');
title('ροπή της θ3');
xlabel('Χρόνος(sec)');
ylabel('Ροπή θ3');
figure;

plot(time, torque_d4(1, :), '*-g');
title('ροπή της d4');
xlabel('Χρόνος(sec)');
ylabel('Ροπή d4');
