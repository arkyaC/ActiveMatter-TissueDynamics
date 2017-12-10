for k=1:10:Nsteps
theta = y(k,2*49+1:3*49);
velX = cos(theta);
velY = sin(theta);
posX=y(k,1:49);
posY = y(k,50:2*49);
plot(posX,posY,'bo','MarkerSize',8,'MarkerFaceColor','b');
axis([0 L 0 L]);
axis square;
grid on;
hold on;
quiver(posX,posY,velX,velY,.65);
pause(0.005);
hold off;

end