% Created by WANG Xuehua 2022/06/29
clc;clear
disp('黑潮-亲潮交汇区涡旋演变过程分析')

%% 载入涡旋中心点、边界信息和生命周期
load('..\data\processeddata\20180701_1231\Tracks\eddy_tracks.mat');
tracks2018.all.day = {tracks.day};
tracks2018.all.lat = {tracks.lat};
tracks2018.all.lon =  {tracks.lon};
tracks2018.all.shapes = {tracks.shapes};
tracks2018.all.type = {tracks.type};
tracks2018.all.atype = cellfirst({tracks.type});
tracks2018.all.lifetime = cellminus({tracks.day});
tracks2018.all.alat = cellfirst({tracks.lat}); % 大致范围
tracks2018.all.alon = cellfirst({tracks.lon});

% 筛选黑潮-亲潮交汇区数据 （145-165E, 38-50N）
flag = find(tracks2018.all.alat>=38&tracks2018.all.alat<=50&tracks2018.all.alon>=145&tracks2018.all.alon<=165);%筛选范围
tracks2018.all.atype=tracks2018.all.atype(flag); % 气旋类型
tracks2018.all.lifetime=tracks2018.all.lifetime(flag); % 气旋生命周期
tracks2018.all.day=tracks2018.all.day(flag); % 存活时段
tracks2018.all.shapes=tracks2018.all.shapes(flag); % 边界信息
tracks2018.all.lat=tracks2018.all.lat(flag); % 中心
tracks2018.all.lon=tracks2018.all.lon(flag);
tracks2018.all.type = tracks2018.all.type(flag);

% 统计研究区域内气旋和反气旋涡数据
CEnum = find(tracks2018.all.atype==1); % 气旋数量
AEnum = find(tracks2018.all.atype==-1); % 反气旋数量
CElt = tracks2018.all.lifetime(find(tracks2018.all.atype==1)); % 气旋生命周期 lifetime
AElt =  tracks2018.all.lifetime(find(tracks2018.all.atype==-1)); % 反气旋生命周期

% % 直方统计图(fig1)
% set(gcf,'Position',[100 50 700 250]);
% hold on;
% edges=[1:1:50]
% h1 = histogram(CElt,edges,'facecolor',[0 0 0]);
% h2 = histogram(AElt,edges,'facecolor',[1 1 1]);
% legend({'CE','AE' });
% box on
% xlabel('生命周期/天','fontweight','bold');
% ylabel('数量','fontweight','bold');
% xlim([0,50])
% set(gca,'xtick',[0:5:50])
% % a=tabulate(tracks2018.all.lifetime(:)); % 统计频率
clear h1 h2 
%% 载入SLA数据
file='..\data\originaldata\sla_20180701_1231_lon_144.5_165.5_lat_37.5_50.5.nc';
ncdisp(file)
SLA2018.all.latitude = ncread(file,'latitude');
SLA2018.all.longitude = ncread(file,'longitude');
SLA2018.all.u = ncread(file,'ugos'); % 流速 m/s
SLA2018.all.v = ncread(file,'vgos');
SLA2018.all.sla = ncread(file,'sla'); % 海表面高度异常值
SLA2018.all.time = ncread(file,'time'); % days since 1950-01-01 00:00:00
SLA2018.all.day = SLA2018.all.time+datenum(1950,1,1) - datenum(2018,7,1)+1;% 匹配老师给的数据日期
clear file

%% 匹配SLA数据和u v数据至涡旋中心点和边界
for i = 1:length(tracks2018.all.day)
    thiseddyday = cell2mat(tracks2018.all.day(i)); % 日期
    thiseddylon = cell2mat(tracks2018.all.lon(i)); % 经度
    thiseddylat = cell2mat(tracks2018.all.lat(i)); % 纬度
    thiseddyedge = tracks2018.all.shapes(i); % 边界信息
    t1=[];
    t2=[];
    t3=[];
    t4=[];
    distance_ab=[];
    for j = 1:length(thiseddyday)
        d = thiseddyday(j);
        a1 = thiseddylon(j);
        a2 = thiseddylat(j);
        b1 = thiseddyedge{1}{j}(1,:); % 经度
        b2 = thiseddyedge{1}{j}(2,:); % 纬度
        distance_ab(j) = mean(distance(a1,a2,b1,b2)/180.*pi.*6371);% 涡旋半径的定义是涡旋边界上各点到涡心的距离平均值。
        [lon2,lat2] = meshgrid(SLA2018.all.longitude,SLA2018.all.latitude);
        eu = interp2(lon2,lat2,SLA2018.all.u(:,:,d)',b1,b2,'spline'); % 匹配的边界流速
        ev = interp2(lon2,lat2,SLA2018.all.v(:,:,d)',b1,b2,'spline');
        csla =interp2(lon2,lat2,SLA2018.all.sla(:,:,d)',a1,a2,'spline'); % 匹配的中心SLA
        esla =interp2(lon2,lat2,SLA2018.all.sla(:,:,d)',b1,b2,'spline'); % 匹配的边界SLA
        t1(j) = csla;
        t2(j) = mean(esla); % 边界sla均值
        t3(j) = mean(esla)-csla; % 振幅 SLA偏差
        % 计算相对涡度
        [lonm,latm]=meshgrid(b1,b2);
        [um,vm]=meshgrid(eu,ev);
        curlz=curlz_atmos(lonm,latm,um,vm);
        t4(j)=mean(curlz(:),'omitnan'); % 边界相对涡度均值
        
    end
    tracks2018.all.radius{i} = distance_ab; % 半径
    tracks2018.all.centersla{i} = t1;
    tracks2018.all.edgesla{i} = t2;
    tracks2018.all.Resla{i} = t3; % 振幅 小于0表示反气旋，下降流，暖流，红色 AE；大于0表示气旋 CE 蓝色
    tracks2018.all.Recurlz{i} = t4;  % 相对涡度
    tracks2018.all.mappedday{i} = mapminmax(cell2mat(tracks2018.all.day(i))', 0, 1); % 标准化的生命周期
end
clear t1 t2 t3 t4 a1 a2 AElt CElt b1 b2 d distance_ab 
clear thiseddyday thiseddyedge thiseddylat thiseddylon i j
clear um vm  eu ev flag lon2 lat2 latm lonm curlz esla csla

%% 涡旋性质的演变 （半径、振幅、相对涡度在标准化生命周期内的平均演变情况）
%最终生成个矩阵
ceer=[]; % eddy radius
ceec=[]; % eddy curl
cees=[]; % eddy relatively SLA
ceed=[]; % 标准化的生命周期
aeer=[]; % 反 radius
aeec=[]; % 反 curl
aees=[]; % 反 relatively SLA
aeed=[]; % 反 标准化的生命周期
scday=[];saday=[];
sclat=[];sclon=[];salat=[];salon=[];
for i = 1:length(tracks2018.all.day)
    if tracks2018.all.atype(i)==1
        ceed =[ceed,cell2mat(tracks2018.all.mappedday(i))];% 标准化的生命周期
        ceer =[ceer,cell2mat(tracks2018.all.radius(i))]; % 气旋半径
        cees =[cees,cell2mat(tracks2018.all.Resla(i))]; % 气旋SLA偏差
        ceec =[ceec,cell2mat(tracks2018.all.Recurlz(i))]; %气旋相对涡度
        sclat=[sclat,cell2mat(tracks2018.all.lat(i))'];
        sclon=[sclon,cell2mat(tracks2018.all.lon(i))'];
        scday=[scday,cell2mat(tracks2018.all.day(i))'];
    elseif  tracks2018.all.atype(i)==-1
        aeed =[aeed,cell2mat(tracks2018.all.mappedday(i))];
        aeer =[aeer,cell2mat(tracks2018.all.radius(i))]; % 气旋半径
        aees =[aees,cell2mat(tracks2018.all.Resla(i))]; % 气旋SLA偏差
        aeec =[aeec,cell2mat(tracks2018.all.Recurlz(i))]; %气旋相对涡度
        salat=[salat,cell2mat(tracks2018.all.lat(i))'];
        salon=[salon,cell2mat(tracks2018.all.lon(i))'];
        saday=[saday,cell2mat(tracks2018.all.day(i))'];
    end
end

CEpos=[ceed',ceer',cees',ceec',sclat',sclon']; % 标准化的生命周期 半径 振幅 相对涡度 纬度 经度
AEpos=[aeed',aeer',aees',aeec',salat',salon'];

%  按标准化的生命周期的大小重新排序，相同周期的特征量取平均值
[bc mc nc]=unique(CEpos(:,1));
for ii=1:length(mc)
    bc(ii,2)=mean(CEpos(nc==ii,2)); % 半径
    bc(ii,3)=mean(CEpos(nc==ii,3)); % 振幅
    bc(ii,4)=mean(CEpos(nc==ii,4)); % 相对涡度
end

[ba ma na]=unique(AEpos(:,1));
for ii=1:length(ma)
    ba(ii,2)=mean(AEpos(na==ii,2));
    ba(ii,3)=mean(AEpos(na==ii,3));
    ba(ii,4)=mean(AEpos(na==ii,4));
end
mapedday2=linspace(0,1,20); % 将各个涡的周期在[0, 1]区间内平均分成20个时间段
% 将各涡旋的特征量插值到这20个时间段内
cenewradius =abs(interp1(bc(:,1),bc(:,2),mapedday2));
aenewradius =abs(interp1(ba(:,1),ba(:,2),mapedday2));
cenewResla =abs(interp1(bc(:,1),bc(:,3),mapedday2));
aenewResla =abs(interp1(ba(:,1),ba(:,3),mapedday2));
cenewRecurlz =abs(interp1(bc(:,1),bc(:,4),mapedday2));
aenewRecurlz =abs(interp1(ba(:,1),ba(:,4),mapedday2));

clear aeec aeed aeer aees bc ba ceec ceed ceed ceer cees salat salon
clear i ii j mc ma nc na sclat sclon saday scday CEnum AEnum

% 画图 (fig2)
set(gcf,'Position',[100 20 1000 700]);
ylabletxt={'半径/km','\deltaSLA/m','相对涡度'};
texttxt={'(a)','(b)','(c)'};
for i=1:3
    subplot(3,1,i)
    y1=[cenewradius;cenewResla;cenewRecurlz];
    scatter(mapedday2,y1(i,:),'filled','b','HandleVisibility','off');
    hold on
    k=polyfit(mapedday2,y1(i,:),5);
    y3=polyval(k,mapedday2);
    plot(mapedday2,y3,'color','b','linewidth',2);
    hold on
    y2=[aenewradius;aenewResla;aenewRecurlz];
    scatter(mapedday2,y2(i,:),'filled','r','HandleVisibility','off');
    hold on
    k=polyfit(mapedday2,y2(i,:),5)
    y3=polyval(k,mapedday2);
    plot(mapedday2,y3,'color','r','linewidth',2);
    if i==3
        ylim([0*10^-9 7*10^-9])
    end
    box on
    set(gca,'Xgrid','on');
    xlabel('标准化的生命周期','fontweight','bold');
    ylabel(ylabletxt{i},'fontweight','bold');
    xli=xlim;
    yli=ylim;
    text(xli(1)+0.03,yli(1)+0.9*(yli(2)-yli(1)),texttxt{i})
end
legend({'CE','AE' });
clear cenewResla aenewResla cenewRecurlz aenewRecurlz aenewradius cenewradius 
clear i k mapedday2 xli yli y1 y2 y3 ylabletxt texttxt

%% 画涡旋半径和振幅的地理分布图(fig3)
% set(gcf,'Position',[100 50 900 500]);
% set(gcf,'color','white');
% for i=1:6
%     % CEpos=[ceed',ceer',cees',ceec',sclat',sclon']; % 标准化的生命周期 半径 振幅 相对涡度 纬度 经度
%     % AEpos=[aeed',aeer',aees',aeec',salat',salon'];
%     if i==1|i==3|i==5
%         A = CEpos;
%     else
%         A = AEpos;
%     end
%     subplot(3,2,i)
%     [C,ia,ic] = unique(A(:,5:6), 'rows');
%     d = A(ia,:);
%     x=d(:,6);
%     y=d(:,5);
%     if i==1|i==2
%         v=abs(d(:,2)); %半径
%     elseif i==3|i==4
%         v=abs(d(:,3)); %振幅
%     else
%         v=abs(d(:,4)); %相对涡度
%     end
%     
%     [xq,yq] = meshgrid(140:1:170, 36:1:55); %将研究海区划分成1°×1°网格
%     vq = griddata(x,y,v,xq,yq); % , 并计算各网格内所有移动经过的涡旋的平均半径、振幅、相对涡度
%     m_proj('equidist', 'lon',[146 164],'lat',[39 49], 'aspect', 0.5);
%     m_pcolor(xq,yq,vq);
%     m_coast('patch', [.7 .7 .7])
%     m_coast('linewidth', 1, 'color', 'b');%画出海岸线
%     m_grid('box','fancy');%添加格网;
%     colorbar
%     colormap(jet);
%     hold on
%     xlabel('经度','FontWeight','bold')
%     ylabel('纬度','FontWeight','bold')
%     % caxis([0 0.4])
%     h = colorbar;
% end
clear A ans C d h i ia ic v vq x xq y yq

%% 画半径概率密度图 (fig4)
% figure
% set(gcf,'Position',[100 50 700 250]);
% % CEpos=[ceed',ceer',cees',ceec',sclat',sclon']; % 标准化的生命周期 半径 振幅 相对涡度 纬度 经度
% % AEpos=[aeed',aeer',aees',aeec',salat',salon'];
% %气旋
% y1 = CEpos(:,2);
% y1min=min(CEpos(:,2));
% y1max=max(CEpos(:,2));
% %反气旋
% y2 = AEpos(:,2);
% y2min=min(AEpos(:,2));
% y2max=max(AEpos(:,2));
% x=linspace(y2min,y2max,20); % 将最大最小区间分成20个等分点(19等分),
% yy1=hist(y1,x); % 计算各个区间的个数
% yy1=yy1/length(y1); % 计算各个区间的个数
% yy2=hist(y2,x);
% yy2=yy2/length(y2);
% % bar(x,yy);%画出概率密度分布图
% scatter(x,yy1,25,'b','LineWidth',1.5);
% hold on;
% scatter(x,yy2,25,'r','LineWidth',1.5);
% xlabel('半径/km','fontweight','bold');
% ylabel('概率密度','fontweight','bold');
% set(gca,'ygrid','on');
% box on
% legend({'CE','AE' });
clear x y1 y1max y1min y2 y2max y2min yy1 yy2

%% 涡旋的振幅、相对涡度与半径的关系 
% 按照半径排序
[bc mc nc]=unique(CEpos(:,2));
for ii=1:length(mc)
bc(ii,2)=mean(CEpos(nc==ii,3));
bc(ii,3)=mean(CEpos(nc==ii,4));
end

[ba ma na]=unique(AEpos(:,2));
for ii=1:length(ma)
ba(ii,2)=mean(AEpos(na==ii,3));
ba(ii,3)=mean(AEpos(na==ii,4));
end
mapedrad2=linspace(10,110,20);%将半径从小到大间隔5km划分成若干个区间,
cenewResla =abs(interp1(bc(:,1),bc(:,2),mapedrad2));
aenewResla =abs(interp1(ba(:,1),ba(:,2),mapedrad2));
cenewRecurlz =abs(interp1(bc(:,1),bc(:,3),mapedrad2));
aenewRecurlz =abs(interp1(ba(:,1),ba(:,3),mapedrad2));

% % 画图(fig5)
% figure
% set(gcf,'Position',[100 50 700 400]);
% ylabletxt={'\deltaSLA/m','相对涡度'};
% texttxt={'(a)','(b)'};
% for  i=1:2
% subplot(2,1,i)
% y1=[cenewResla;cenewRecurlz];
% scatter(mapedrad2,y1(i,:),'filled','b','HandleVisibility','off');
% hold on
% k=polyfit(mapedrad2,y1(i,:),5)
% y3=polyval(k,mapedrad2);
% plot(mapedrad2,y3,'color','b','linewidth',2);
% hold on
% y2=[aenewResla;aenewRecurlz];
% scatter(mapedrad2,y2(i,:),'filled','r','HandleVisibility','off');
% hold on
% k=polyfit(mapedrad2,y2(i,:),5)
% y3=polyval(k,mapedrad2);
% plot(mapedrad2,y3,'color','r','linewidth',2);
% box on
% set(gca,'Xgrid','on');
% xlabel('半径/km','fontweight','bold');
% ylabel(ylabletxt{i},'fontweight','bold');
% xli=xlim;
% yli=ylim;
% text(xli(1)+0.2,yli(1)+0.9*(yli(2)-yli(1)),texttxt{i})
% end
% legend({'CE','AE' });
clear cenewResla aenewResla cenewRecurlz aenewRecurlz aenewradius cenewradius 
clear ba bc i ii k ma mc nc texttxt xli y1 y2 y3 ylabletxt yli na mapedrad2


%%  研究区域内各涡海平面动力高度异常(fig6)
% j=0
% figure('position',[333,44,600,300]);
% for i=1:length(tracks2018.all.lifetime)
%     lt=tracks2018.all.lifetime(i)
%     if lt>=56&lt<63 % 周
%         j=j+1
%         %  subplot(7,2,1)
%         x0  = cell2mat(tracks2018.all.day(i))';
%         x=x0-x0(1);
%         y=cell2mat(tracks2018.all.Resla(i));
%         plot(x,y,'linewidth',1,'color','k');
%         hold on
%     end
%     set(gca,'linewidth',1,'fontsize',10,'fontweight','bold');
%     set(gca,'XTick',[0:7:63]);
%     %     set(gca,'XTickLabel',{'0','7','14','21','});
%     %     set(gca,'YTick',[-9*10^(-8):3*10^(-8):9*10^(-8)]);
%     set(gca,'XLim',[0 63]);
%     set(gca,'YLim',[-0.2 0.2]);
%     xlabel('Day(8-9 weeks)','fontsize',12,'fontweight','bold');
%     ylabel('\deltaSLA/m','fontsize',12,'fontweight','bold');
%     grid on
% end
clear i j lt x x0 y