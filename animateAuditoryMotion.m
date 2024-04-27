function [yy,p] = animateAuditoryMotion(p,s,x,y,n)

if ~exist('n','var')
   n=1;
end



yy = makeAuditoryMotion(p,s,x,y);

figure(1)
clf


hold on

%'Head'
plot(0,0,'ko','MarkerSize',10,'MarkerFaceColor','y');
axis equal

plot(x(1:200:end),y(1:200:end),'r:');
h = plot(x(1),y(1),'ko','MarkerFaceColor','r','MarkerSize',5);
xlabel('X (meters)');
ylabel('Y (meters)');
yy = yy/max(yy(:));


%loop n times

for i=1:n
    %play the sound
    sound(yy,p.Fs);
    %animate the object
    tic
    while toc<p.dur
        %Find index that matches closest to the current time
        id = ceil(toc*p.Fs);
        set(h,'XData',x(id));
        set(h,'Ydata',y(id));
        drawnow
    end
end