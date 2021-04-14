
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Multiple Restart Hill Climbing Algorithm for a bidimentional function for
%a minimization problem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHASE 1- Setting up for the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
func = @(X,Y) 0.5+((sin(sqrt(X.^2+Y.^2))).^2-0.5/1+0.001.*(X.^2+Y.^2).^2);
%Visual representation in 3D: %
scrsz = get(0,'ScreenSize');
f1 = figure(1);
surf(X,Y, 0.5+((sin(sqrt(X.^2+Y.^2))).^2-0.5/1+0.001.*(X.^2+Y.^2).^2))
%labels
xlabel('x1')
ylabel('x2')
zlabel('f')
colormap(bone);
%Implementing graphical limits %
limits = [-2.1 2.1 -2.1 2.1 0 1.5];
axis(limits)
colorbar
hold on %retains the current plot
pause(1) %pauses during 1 second

%Defining bounds of the function, its maximum radius, and its moving point
%neighborhood radius
VarLBounds=[-2.048 -2.048];
VarUBounds=[+2.048 +2.048];
max_radius = VarUBounds - VarLBounds;
radius = 100;
%Initializing first x and y points in the function! %
xy_current = VarLBounds(1,:)+ rand(1,2) .* max_radius(1,:);
% Create current solution
s_current = func(xy_current(1),xy_current(2));
%Defining a color vector to operate on, so different colors can be shown
ColorMarker = [1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHASE 2- Running MRHC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ready to start!
disp('--------Multiple Restart Hill Climbing----------')
disp('Hill Climbing in progress')
%Plot the first obtained x and y points from line 49%
plot3(xy_current(1),xy_current(2),s_current,'ok','MarkerSize',5,'MarkerFaceColor',ColorMarker);
%Saving the points so we can later build a graphical representation %
x_plot(1) = xy_current(1);
y_plot(1) = xy_current(2);
f_plot(1) = s_current;
res = 1; %iterator to track renicialization
%Defining the iterator that tracks repetition (function outputs same values)%
%(stagnation criteria) %
it_rep = 0;
res_point(1) = 0; %Creating a vector to save the iteration where a restart
%is performed
%For Cycle
for j=1:1000
    %New points of x/y in the neighborhood of current points
    xy_new(j,:) = xy_current + (((VarUBounds-VarLBounds)/radius).*randn(1,2));
    %Truncating the values to the limits of the function (for both x and y)
    if xy_new(j,1)<VarLBounds(1,1)
       xy_new(j,1)=VarLBounds(1,1);
    end
    if xy_new(j,2)<VarLBounds(1,2)
       xy_new(j,2)=VarLBounds(1,2);
    end
    if xy_new(j,1)>VarUBounds(1,1)
       xy_new(j,1)=VarUBounds(1,1);
    end
    if xy_new(j,2)>VarUBounds(1,2)
       xy_new(j,2)=VarUBounds(1,2);
    end
    % Create new solution
    s_new = func(xy_new(j,1),xy_new(j,2));
    if (s_new < s_current)%New solution is better than the current?..
        xy_current(1) = xy_new(j,1);
        xy_current(2) = xy_new(j,2);
        s_current = s_new; %Updates solution
        %Plot the new points
        plot3(xy_current(1),xy_current(2),s_current,'ok','MarkerSize',5,'MarkerFaceColor',ColorMarker);
        pause(0.1) %Pause introduces a delay before new points are plotted
        it_rep=0; %Stagnation criteria is set to 0 since the algorithm successfully found new points
    end
    % Checking if the new solution is equal to the old one, so it_rep can
    %keep track if the hill climbing is stuck in a local minimum
    if(abs(s_new - s_current) < 0.0001 && it_rep ~= 5)
        it_rep = it_rep + 1;
    end
    %Renitialization %
    %if the solution is the same about 5 more cycles, than x and y are
    %randomly generated (like in the beginning)
    if(it_rep == 5)
        res = res + 1;
        res_point(res) = j; %res_points is a vector that saves the j
                          %(iterator) points in which there was a restart
                          %for plotting purposes
        xy_current = VarLBounds(1,:)+ rand(1,2) .* max_radius(1,:);
        resx_plot(res) = xy_current(1);
        resy_plot(res) = xy_current(2);
        s_current = func(xy_current(1),xy_current(2));
        sol_plot(res) = s_current;
        ColorMarker(1) = rand;
        ColorMarker(2) = rand;
        ColorMarker(3) = rand;
        plot3(xy_current(1),xy_current(2),s_current,'d','MarkerSize',5,'MarkerFaceColor',ColorMarker);
        it_rep=0;
    end
    %Saving all the generated points
    x_plot(j) = xy_current(1); 
    y_plot(j) = xy_current(2);
    f_plot(j) = s_current;
    %increment iterator
    j=j+1;
end
%Set the graphics free
hold off
%Graphical representation of the values of x, y and the solution during the
%the iterations in the cycle
%For the algorithm to be successfull, all values should converge to the
%point of the coordinates x,y (0,0) where the value of the function
%(solution) is 0.
f2 = figure('Position',[25 380 670 220]); % [left, bottom, width, height]
subplot(2,1,1)
plot(x_plot,'r-','LineWidth',1)
xlabel('iterations')
ylabel('current x')
yline(0,'--k','LineWidth',1.1);
subplot(2,1,2)
plot(y_plot,'r-','LineWidth',1)
xlabel('iterations')
ylabel('current y')
yline(0,'--k','LineWidth',1.1);
f3 = figure('Position',[25 110 670 180]);
plot(f_plot,'r-','LineWidth',1)
xlabel('iterations')
ylabel('current solution')
yline(0,'--k','LineWidth',1.1);
%in this last plot, there will be representations of the restart points
%if any occured
for cont = 2: length(res_point)
    lines = xline(res_point(cont),'--b','LineWidth',1.3,'DisplayName','Restart points');
end
if res >= 2
    legend(lines,'Restart points');
    f4 = figure('Position',[700 135 540 460]);
    ezcontourf(func,[-2.048 2.048 -2.048 2.048])
    colormap(bone);
    hold on
    plot(resx_plot(2:end),resy_plot(2:end),'d','MarkerSize',15,'MarkerFaceColor','b','MarkerEdgeColor','white')
    title('Restart points')
    hold off
end
