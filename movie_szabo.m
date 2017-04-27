%fid = fopen('noswitch_horseshoe_18_05.dat');

fid = fopen('movie_dump.txt');

[data] = textscan(fid,'%f %f %f %f %f %f');

fclose(fid);

 

%Initializing data

 

time = data{:,1};

time = int32(time);

r(:,1) = data{:,3};

r(:,2) = data{:,4};

v(:,1) = data{:,5};

v(:,2) = data{:,6};

tmax = max(time)-1;

tmin = min(time);

 

 

Nc = length(time(time==tmax));

 

xmat = NaN(tmax,Nc,2);

xvec = NaN(tmax,Nc,2);

for i=1:tmax

    n = length(r(time==i,1));

    xmat(i:i,1:n,1)=r(time==i,1);

    xmat(i:i,1:n,2)=r(time==i,2);

    xvec(i:i,1:n,1)=v(time==i,1);

    xvec(i:i,1:n,2)=v(time==i,2);

end

 

 

figure

%axis tight

%set(gca,'nextplot','replacechildren');

 

%vidObj = VideoWriter('T_adh.avi');

%vidObj.FrameRate = 20;

%open(vidObj);

 

% Record the movie

n=0;

qscale=0.7;

for j = tmin:1:tmax

%for j = 1:1:tmax-1

    

    n=n+1;

    

    plot(xmat(j,:,1),xmat(j,:,2),'ob','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','b');

    hold on;

    h= quiver(xmat(j,:,1),xmat(j,:,2),xvec(j,:,1),xvec(j,:,2),0,'k','LineWidth',1);

    

    hU = get(h,'UData') ;

    hV = get(h,'VData') ;

    set(h,'UData',qscale*hU,'VData',qscale*hV)

    

    axis equal

    xlim([0 20]);

    ylim([0 20]);

    set(gca,'XTick',[])

    set(gca,'YTick',[])

    hold off;

    

    F(:,n:n) = getframe;

    

    %writeVideo(vidObj,F);

    

end

%close(vidObj);

movie(F,100);

