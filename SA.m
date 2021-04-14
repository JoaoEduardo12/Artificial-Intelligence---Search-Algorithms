%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulated Annealing Algorithm for a bidimentional function for
%a minimization problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHASE 1 -Setting up for the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note the 1st phase is equal to the 1st phase to the MRHC with slight
%changes on lines 62

%Close all the windows, clear working space, shut off warnings
clc
clear
close all 
warning('off','all');
% Function Initialization and Representation %
X = -2.048:0.08:2.048; %vector from -2.048 to 2.048 with increments of 0.08
Y = X; %Y as the same vector
[X,Y] = meshgrid(X); %returns square grid coordinates with the size of X
%Defining the function func: %
func = @(X,Y) 0.5+((sin(sqrt(X.^2+Y.^2))).^2-0.5/1+0.001*(X.^2+Y.^2).^2);
%Visual representation in 3D: %
f1 = figure(1);
surf(X,Y, 0.5+((sin(sqrt(X.^2+Y.^2))).^2-0.5/1+0.001*(X.^2+Y.^2).^2));
%labels
xlabel('x1')
ylabel('x2')
zlabel('f')
%Implementing graphical limits %
limits = [-2.1 2.1 -2.1 2.1 0 1.5];
axis(limits)
%Colorbar and retain current plot %
colormap(bone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2 = figure(2);
% ezcontourf(func,[-2.048 2.048 -2.048 2.048]);
% colormap
% hold on
% plot(0,0,'p','MarkerSize',15,'MarkerFaceColor','yellow','MarkerEdgeColor','k')
% plot(-2.048,-2.048,'p','MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','k')
% plot(2.048,2.048,'p','MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','k')
% plot(2.048,-2.048,'p','MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','k')
% plot(-2.048,2.048,'p','MarkerSize',15,'MarkerFaceColor','r','MarkerEdgeColor','k')
colorbar %can add different color palettes
hold on %retains the current plot
pause(1) %pauses during 1 second
%Defining function's boundaries
VarLBounds=[-2.048 -2.048];
VarUBounds=[+2.048 +2.048];
max_radius = VarUBounds - VarLBounds;
radius = 30;
%Initializing first x and y points in the function! %
xy_current = VarLBounds(1,:)+ rand(1,2) .* max_radius(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHASE 2- Running SA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ready to start!
disp('--------------Simulated Annealing---------------')
%Note: Our base variables are set, however we have to set SA's specific exclusive
%variables
ColorMarker = [1 0 0];
% Set initial temperature
t_initial = 4;
% Set final temperature
t_final = 0.3;
% Create current solution
s_current = func(xy_current(1),xy_current(2));
%Plot the first point
plot3(xy_current(1),xy_current(2),s_current,'ok','MarkerSize',4,'MarkerFaceColor',ColorMarker);
%Saving the first point
x_plot(1) = xy_current(1);
y_plot(1) = xy_current(2);
f_plot(1) = s_current;
% Set the temperature
t = t_initial;
%Saving the first value of temperature
t_plot(1) = t;
%Defining the probability funtion here so we can save the first value
%(this is again defined within the algorithm so it can change its value
%during the run)
p = 0.8;
%Save the first point
p_plot(1) = p;
disp(['Starting temperature is ',num2str(t_initial)])
disp('Simulation starting...')
%We do need an extra iterator to keep track of the total ammount of
%iterations, so we can save the values of x, y, s, temperature and
%probability
i = 2;
delta_values(1) = 1;
% Temperature Cycle
while t > t_final %as long as the temperature is bigger than the final
                  %temperature, the SA keeps running
    iter = 1;     %Setting the iterations
    % Metropolis Cycle
    while iter < 100  %As long as the iterations dont surpass 200, the
                      %inner cycle keeps running
        %Generate new points of x and y from the neighborhood of the
        %current points
        xy_new(iter,:) = xy_current + (((VarUBounds-VarLBounds)/radius).*randn(1,2));
        %Truncating the values (x,y) to the limits of the function   
        if xy_new(iter,1)<VarLBounds(1,1)
           xy_new(iter,1)=VarLBounds(1,1);
        end
        if xy_new(iter,2)<VarLBounds(1,2)
           xy_new(iter,2)=VarLBounds(1,2);
        end
        if xy_new(iter,1)>VarUBounds(1,1)
           xy_new(iter,1)=VarUBounds(1,1);
        end
        if xy_new(iter,2)>VarUBounds(1,2)
           xy_new(iter,2)=VarUBounds(1,2);
        end    
        % Create new solution
        s_new = func(xy_new(iter,1),xy_new(iter,2));
        % Calculate the probability p with s_new and s_current using a
        % boltzmann distribution, the value t will continue to change in
        % the outer cycle (Temperature Cycle), so that p gets smaller and
        % smaller.
        delta_s = s_new - s_current;
        delta_values(iter) = delta_s;
        delta_s_avg = (sum(abs(delta_values),[2:length(delta_values)])/length(delta_values)-1);
        p = exp(-abs(delta_s_avg)/t); %probability function based on boltzmann distribution
        % Comparing the values of both new and old solutions
        if s_new < s_current %New solution is better than the current?..
            s_current = s_new;%.. than the new solution becomes the current
            xy_current(1) = xy_new(iter,1);
            xy_current(2) = xy_new(iter,2);
        % Another chance for s_new to be accepted even though its a worse
        % solution for minimization (since it has a bigger value and did
        % not run in the if statement line 151).
        elseif rand < p   %a randomly generated number between 0 and 1 is
                          %smaller than p?..
            s_current = s_new; %.. then we accept the new solution
            xy_current(1) = xy_new(iter,1);
            xy_current(2) = xy_new(iter,2);
        end
        plot3(xy_current(1),xy_current(2),s_current,'ok','MarkerSize',4,'MarkerFaceColor',ColorMarker);
        pause(0.1)
        %Save all the points x,y, the solutions and the temperature and
        %the probability function
        x_plot(i) = xy_current(1);
        y_plot(i) = xy_current(2);
        f_plot(i) = s_current;
        t_plot(i) = t;
        p_plot(i) = p*100;
        i = i + 1;
        iter = iter + 1;
    end
    % Cooling function for temperature decrease
    t = 0.87 * t;
    if (p * 100) < 80
        ColorMarker = [0.6 0.3 0.3];
    end
    if (p * 100) < 55
        ColorMarker = [0.4 0.5 0.5];
    end
    if (p * 100) < 30
        ColorMarker = [0.2 0.75 0.75];
    end
    if (p * 100) < 10
        ColorMarker = [0 1 1];
    end
    disp(['Current temperature is now: ',num2str(t)])
    disp(['Current probability of accepting worse solutions is now: ',num2str(p * 100),'%'])
end
hold off %Sets the graphs free
%Graphical representation of the values of x, y and the solution during the
%the iterations in the cycle
%For the algorithm to be successfull, all values should converge to the
%point of the coordinates x,y (0,0) where the value of the function
%(solution) is 0.
figure('Position',[100 120 500 440])
subplot(3,1,1)
plot(x_plot,'b-','LineWidth',1)
yline(0,'--k','LineWidth',1.1);
xlabel('current x')
subplot(3,1,2)
plot(y_plot,'b-','LineWidth',1)
yline(0,'--k','LineWidth',1.1);
xlabel('current y')
subplot(3,1,3)
plot(f_plot,'b-','LineWidth',1)
yline(0,'--k','LineWidth',1.1);
xlabel('current solution value')

figure('Position',[700 120 500 440])
subplot(2,1,1)
plot(t_plot,'b-','LineWidth',1.1)
xlabel('iterations')
ylabel('Temperature')
subplot(2,1,2)
plot(p_plot,'b-','LineWidth',1.1)
xlabel('iterations')
ylabel('Probability')
