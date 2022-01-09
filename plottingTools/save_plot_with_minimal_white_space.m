%% Save plot with minimal white space without removing part of the figure
% https://www.mathworks.com/matlabcentral/answers/466557-save-plot-with-minimal-white-space-without-removing-part-of-the-figure


% Set up figure size before moving forward.
% Figure units must be pixels (default).
fh = figure('Units','Pixels'); 
% Set axes is fill figure
ax = axes('position',[0 0 1 1]);
% create sample data 
x = -10:10; 
y = -10:10; 
[X,Y] = meshgrid(x,y); 
Z = X.^2 + Y.^2; 
% plot 
contourf(ax,X,Y,Z); % Use axis handle!
cbh = colorbar(ax); % Use axis handle!
cbh.Position([2,4]) = [.015,.97];  % decrease vertical size of colorbar
set(ax,'xtick',[])
set(ax,'ytick',[])
axis(ax,'equal')
% Set colorbar units and axis units to match figure units
cbh.Units = fh.Units; 
ax.Units = fh.Units; 
% Move axis all the way to the left
ax.Position(1) = -(ax.Position(3)-ax.Position(4))/2 + 1;
% Move colorbar to the left (the +5 at the end is pixel space between plot and colorbar)
cbh.Position(1) = sum(ax.Position([1,3])) + ax.Position(1) + 5;
% Now shrink the width of the figure; the *1.5 at the end is to factor in the y tick marks on the colorbar. 
fh.Position(3) = sum(cbh.Position([1,3]))+cbh.Position(3)*1.5; 
% Turn off figure border (if that's what you want to do)
fh.MenuBar = 'none';
fh.ToolBar = 'none';

% Comment: Any chance you can point me to a resource or explain what the
% position values mean for the axis, figure, and color handles? Just so I
% can edit this later if need be.  I assume figure handle is the actual
% outer figure window, axis is specifically the plot, and of course the
% colorbar handle is for the colorbar? Is the Position vector organized as
% [xmin ymin xmax ymax]?
% Comment: Sure!  Check out the properties listed in the "position" group
%          for axes: 
%            
%          https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html#d117e68467
%          The tricky part that isn't so clear in the documentation is that
%          once you make the axes equal, the visible border of the axes
%          doesn't really indicate its true position.  For example, in the
%          figure I posted, the bottom, left corner is not at position
%          (0,0) in the figure.  It's lateral position is actually
%          negative.  Even though the axis appears square it's actual
%          invisible width is wider than its height.  To compensate for
%          that, I moved the axes leftward by 1/2 of the difference between
%          the width and height; hence; 
%          ax.Position(1) = -(ax.Position(3)-ax.Position(4))/2 + 1;
%
%         Also recall that we're workng in pixel units rather than
%         normalized units.
%
%%   More Answers (1)

% If you are not trying to develop a flexible application and just need a
% quick solution, just adjust the numbers in the following:

figure()
set(gca,'position',[0.07 0.07 0.92 0.88]);

% done.


