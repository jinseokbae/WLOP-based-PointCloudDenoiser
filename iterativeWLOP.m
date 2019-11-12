close all
clear all

%%Used a segmented view of CMU dataset of 160906 pizza1 
%% Load Data
ptCloud_ori = pcread('segmented1.ply');
gridStep = 3;
ptCloud = pcdownsample(ptCloud_ori,'gridAverage',gridStep);

m = ptCloud.Count;
pointSet = ptCloud.Location;
%% Random Point Set
XLimits = ptCloud.XLimits;
YLimits = ptCloud.YLimits;
ZLimits = ptCloud.ZLimits;

d_bb = ((XLimits(2)-XLimits(1))^2+(YLimits(2)-YLimits(1))^2+...
    (ZLimits(2)-ZLimits(1))^2)^0.5;
h = 4*(d_bb)^0.5;
n = 1024; % Number of arbitrary points at the initial point

randomXset = rand(n,1)*(XLimits(2)-XLimits(1)) + XLimits(1);
randomYset = rand(n,1)*(YLimits(2)-YLimits(1)) + YLimits(1);
randomZset = rand(n,1)*(ZLimits(2)-ZLimits(1)) + ZLimits(1);

randomSet = [randomXset randomYset randomZset];
%randomSet = pointSet(randi(m,n,1),:); 


%% iteration

iterNum = 25;

R_v = repmat(pointSet,1,m)-repmat(reshape(pointSet',1,[]),m,1);
realR_v = zeros(m,m);
for i = 1:m
    realR_v(:,i) = vecnorm(R_v(:,3*i-2:3*i),2,2);
end
R_v = [];
v = sum(theta(realR_v,h),2);
realR_v = [];

for k = 1:iterNum    
    newSet = zeros(size(randomSet));
    
    R_w = repmat(randomSet,1,n)-repmat(reshape(randomSet',1,[]),n,1);
    realR_w = zeros(n,n);
    for i = 1:n
        realR_w(:,i) = vecnorm(R_w(:,3*i-2:3*i),2,2);
    end
    w = sum(theta(realR_w,h),2);
    realR_w = [];
    
    alpha = zeros(n,m);
    beta = zeros(n,n);
    for i = 1:n
        for j = 1:m
            r = randomSet(i,:)-pointSet(j,:);
            r = norm(r);
            alpha(i,j) = theta(r,h)/(r*v(j));
        end
        
        for j = 1:n
            if(i ~= j)
                r = randomSet(i,:)-randomSet(j,:);
                r = norm(r);
                beta(i,j) = -theta(r,h)*w(j)/r;
            end
        end
    end
    sumalpha = sum(alpha,2);
    sumbeta = sum(beta,2); 
    for i = 1:n
        for j = 1:m
            newSet(i,:) = newSet(i,:) + pointSet(j,:)*alpha(i,j)/sumalpha(i);
        end
        
        for j = 1:n
            if(i~=j)
                newSet(i,:) = newSet(i,:)+0.45*R_w(i,3*j-2:3*j)*beta(i,j)/sumbeta(i);
            end
        end
    end
    
    randomSet = newSet;
    
    
    
    %randomSet = newSet;
    fprintf("iter %d has done\n",k);
    figure(1)
    hold on
    if k==1
        pcshow(ptCloud_ori);
        set(gcf,'color','w')
        axis off
        hold on
    scatter3(newSet(:,1)+50*k,newSet(:,2),newSet(:,3),'.','MarkerFaceColor',[0 .75 .75])
    set(gcf,'color','w')
    axis off
    elseif k<4
    hold on
    scatter3(newSet(:,1)+50*k,newSet(:,2),newSet(:,3),'.','MarkerFaceColor',[0 .75 .75])
    set(gcf,'color','w')
    axis off
    elseif k==20
    hold on
    scatter3(newSet(:,1)+50*4,newSet(:,2),newSet(:,3),'.','MarkerFaceColor',[0 .75 .75])
    set(gcf,'color','w')    
    axis off
    end
    
end

%%
figure;
pcshow(pcdownsample(ptCloud,'gridAverage',20))
set(gcf,'color','w')
axis off
hold on
scatter3(newSet(:,1),newSet(:,2),newSet(:,3),'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0])
set(gcf,'color','w')
hold on
xyz = ptCloud_ori.Location;
ptCloud_zzab = pointCloud(xyz+repmat([50 0 0],length(xyz),1));
ptCloud_zzab.Color = ptCloud_ori.Color;
pcshow(ptCloud_ori)
set(gcf,'color','w')
axis off


%% functionize theta, eta

function ret = theta(r,h)
if(length(r)~=1)
ret = exp(-r.^2/((h/4)^2));

else
  ret = exp(-r^2/((h/4)^2));  
end
end

function ret = dervEta(r)
ret = -1;
% or use
% ret = -1;

end