cvmin=-3.13;
cvmax=3.13;
cvbins=79;
Rmin=5000;

for i=1:cvbins
    for j=1:cvbins
        X(j+(i-1)*cvbins,1)=(cvmax-cvmin)/cvbins*(i-0.5) + cvmin;
        X(j+(i-1)*cvbins,2)=(cvmax-cvmin)/cvbins*(j-0.5) + cvmin;
    end
end

Output_Data = importdata('out (copy).dat');
Output_Data(isnan(Output_Data)) = 0;

N(:,1)=Output_Data(:,3);
N(N<Rmin)=Rmin;
Y(:,1)=Output_Data(:,1)./N(:,1);
Y(:,2)=Output_Data(:,2)./N(:,1);
% for i=1:3
% Y((Y(:,1))>10,2)=Y((Y(:,1))>10,2)./Y((Y(:,1))>10,1)*10;
% Y((Y(:,1))>10,1)=10;
% 
% Y((Y(:,2))>10,1)=Y((Y(:,2))>10,1)./Y((Y(:,2))>10,2)*10;
% Y((Y(:,2))>10,2)=10;
% 
% Y((Y(:,1))<-10,2)=Y((Y(:,1))<-10,2)./Y((Y(:,1))<-10,1)*-10;
% Y((Y(:,1))<-10,1)=-10;
% 
% Y((Y(:,2))<-10,1)=Y((Y(:,2))<-10,1)./Y((Y(:,2))<-10,2)*-10;
% Y((Y(:,2))<-10,2)=-10;
% end

  %Y(X(:,2)==0,2)=0;
  %Y(X(:,2)==0,1)=0;
  %Y(X(:,1)==0,2)=0;
  %Y(X(:,1)==0,1)=0;

figure(1)
q=quiver(X(1:end-cvbins,2),X(1:end-cvbins,1),-Y(1:end-cvbins,2),-Y(1:end-cvbins,1),'k','linewidth',1);
q.AutoScale='on';
q.AutoScaleFactor=3;
axis([cvmin cvmax cvmin cvmax]);

G=zeros(cvbins.^2,1);
for i=1:cvbins
    for j=1:cvbins
        G((i-1)*cvbins+j)= G((i-1)*cvbins+j) + sum(Y(1:j,2));
        for k=1:i
            G((i-1)*cvbins+j)= G((i-1)*cvbins+j) + Y((k-1)*cvbins+j,1);
        end
    end
end

for i=1:cvbins
    for j=1:cvbins
        Ysurf(i,j)=X((i-1)*cvbins+j,1);
        Xsurf(i,j)=-X((i-1)*cvbins+j,2);
        Zsurf(i,j)=G((i-1)*cvbins+j);
    end
end

figure(2)
surf(Xsurf,Ysurf,Zsurf);
axis([cvmin cvmax cvmin cvmax]);

for i=1:cvbins
    for j=1:cvbins
        for k=1:3
            for l=1:3
                F((i-1)*9*cvbins+9*(j-1)+l+k*cvbins,1) = Output_Data((i-1)*cvbins+j,1)/Output_Data((i-1)*cvbins+j,3)*2000;
                F((i-1)*9*cvbins+9*(j-1)+l+k*cvbins,2) = Output_Data((i-1)*cvbins+j,2)/Output_Data((i-1)*cvbins+j,3)*2000;
                F((i-1)*9*cvbins+9*(j-1)+l+k*cvbins,3) = 2000;
            end
        end
    end
end

figure(2)
surf(Xsurf,Ysurf,Zsurf);
axis([cvmin cvmax cvmin cvmax]);


        

