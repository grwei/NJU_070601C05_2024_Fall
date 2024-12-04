%--------------------------------------------------
%Problem Set #5 - Organic waste disposal in a river
%--------------------------------------------------

%Here's some sample code to plot the BO and DO curves through the centerline
%of the source location (select index yp accordingly). You can rename variables,
%revise, and improve this code as you wish.

%At the beginning of your code you can create an avi file. Note that if the file
%is open in a viewer like windows media player, this will give an error. Close
%the movie before rerunning your matlab code.
aviobj1 = VideoWriter('DOsag_curve.avi')

%You can set properties of the video before you open it
aviobj1.FrameRate = 2; %sets the frame rate for the movie to 2 frames/second (adjust as you wish to create a nice, easily understandable movie)

%open the avi object
open(aviobj1);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%You can insert the code below into your time loop, but you may want to add
%a modulo command to only collect a frame every few time steps, e.g. with
% if (mod(m,10)==0) ... where m is your time step index. This may be
%especially worthwhile if you start increasing your resolution and decreasing
%your time step, or for the 3D graphics which take more disk space.

fig = figure(1); %this makes sure to grab the same figure at each time step

%plotyy will plot each variable on a different y axis. do_old is DO and c_old
%is the BOD concentration.
[ax,h1,h2]=plotyy(x,do_old(:,yp),x,c_old(:,yp));

%this speeds up the display of the image
set(fig,'DoubleBuffer','on');

%You can set the axis limits for each y axis separately so that the axis
%does not change with time, even if your concentration range changes with
%time. Choose the tick mark intervals that are best for you. Here, Qload
%is the load strength specified by the user.
set(ax(1),'ylim',[4 8],'ytick',[4 4.5 5 5.5 6 6.5 7 7.5 8])
set(ax(2),'ylim',[0 Qload*1.5], 'ytick',[0 Qload*0.25 Qload*0.5 Qload*0.75 Qload*1 Qload*1.25 Qload*1.5])

%Again, this improves display and capture of the image in the figure window.
set(gca,'nextplot','replace','Visible','on')

%Convert the current time to a string for printing in the title
timestr = sprintf('%5.3f',time(m));
%Display this title - you may want to add your name here
title(['DO sag curve along j = ',num2str(yp),' t = ',timestr,'s']);

%Capture the frame
frame = getframe(gcf);
%Add the frame to your avi movie file
writeVideo(aviobj1,frame);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%At the end of your code (after all time stepping is over), close the movie file
close(aviobj1); %close DO sag movie

%You will create a movie using similar commands for the surface plot (use
%surf command) of dissolved oxygen. There you might want to set 'zlim'
%so that the z axis is fixed in time with this for example:
set(gca,'xlim',[0 Lx],'ylim',[0 Ly],'zlim',[5 8],...
  'nextplot','replace','Visible','on')

