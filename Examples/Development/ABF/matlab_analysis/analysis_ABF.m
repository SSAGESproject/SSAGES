cvmin=2.2;
cvmax=10;
cvbins=120;
Rmin = 200;

Output_Data = importdata('out (copy).dat');
%Output_Data=outcopy;
xf = linspace(cvmin,cvmax,size(Output_Data,1));
n = Output_Data(:,2);
n(n<Rmin) = Rmin;
yf = Output_Data(:,1)./n;

%plot(xpmf,Output_Data);

%for i=2:size(xf,2)-1
 %   yf(i) = Output_Data(i-1)+2*Output_Data(i)+Output_Data(i+1);
%end

figure(1)
plot(xf,yf)
axis([2 10 -5 5]);

xpmf = xf;
ypmf = yf;

for i=2:size(xpmf,2)
    ypmf(i) = ypmf(i-1)+yf(i)*((xf(i)-xf(i-1)));
    ypmfc(i) = ypmf(i) - 2*0.519*log(xpmf(i));
end

figure(2)
%plot(xpmf(5:end),-2*ypmf(5:end)+3.1)
plot(xpmf,-ypmf+70,xpmf,-ypmfc+69);
axis([2 10 -4 2]);
legend('No Entropy Correction','Entropy Correction');

% for i= 1:size(Output_Data,1)-1
%     interp(2*i-1,1) = Output_Data(i,1)/2;
%     interp(2*i,1) = (Output_Data(i,1)+Output_Data(i+1,1))/4;
%     interp(2*i-1,2) = Output_Data(i,2)/2;
%     interp(2*i,2) = (Output_Data(i,2)+Output_Data(i+1,2))/4;
% end
% interp(80,1) = Output_Data(40,1)/2;
% interp(80,2) = Output_Data(40,2)/2;
% interp(:,2) = round(interp(:,2));


