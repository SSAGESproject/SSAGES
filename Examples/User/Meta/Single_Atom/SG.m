function f = SG(x,y,xc,yc,height,sigma)
xdiff = x-xc;
ydiff = y-yc;

for(i=1:length(xdiff))
   if(xdiff(i) > pi())
       xdiff(i) = xdiff(i) - 2*pi();
   end
   if(xdiff(i) < -pi())
       xdiff(i) = xdiff(i) + 2*pi();
   end
   
end

for(i=1:length(ydiff))
   if(ydiff(i) > pi())
       ydiff(i) = ydiff(i) - 2*pi();
   end
   if(ydiff(i) < -pi())
       ydiff(i) = ydiff(i) + 2*pi();
   end
end

f = sum(height*exp(-((xdiff).^2 + (ydiff).^2)./(2*sigma^2)));