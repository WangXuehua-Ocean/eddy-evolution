% Created by WANG Xuehua 2022/06/29
clc;clear
disp('�ڳ�-�׳������������ݱ���̷���')

%% �����������ĵ㡢�߽���Ϣ����������
load('..\data\processeddata\20180701_1231\Tracks\eddy_tracks.mat');
tracks2018.all.day = {tracks.day};
tracks2018.all.lat = {tracks.lat};
tracks2018.all.lon =  {tracks.lon};
tracks2018.all.shapes = {tracks.shapes};
tracks2018.all.type = {tracks.type};
tracks2018.all.atype = cellfirst({tracks.type});
tracks2018.all.lifetime = cellminus({tracks.day});
tracks2018.all.alat = cellfirst({tracks.lat}); % ���·�Χ
tracks2018.all.alon = cellfirst({tracks.lon});

% ɸѡ�ڳ�-�׳����������� ��145-165E, 38-50N��
flag = find(tracks2018.all.alat>=38&tracks2018.all.alat<=50&tracks2018.all.alon>=145&tracks2018.all.alon<=165);%ɸѡ��Χ
tracks2018.all.atype=tracks2018.all.atype(flag); % ��������
tracks2018.all.lifetime=tracks2018.all.lifetime(flag); % ������������
tracks2018.all.day=tracks2018.all.day(flag); % ���ʱ��
tracks2018.all.shapes=tracks2018.all.shapes(flag); % �߽���Ϣ
tracks2018.all.lat=tracks2018.all.lat(flag); % ����
tracks2018.all.lon=tracks2018.all.lon(flag);
tracks2018.all.type = tracks2018.all.type(flag);

% ͳ���о������������ͷ�����������
CEnum = find(tracks2018.all.atype==1); % ��������
AEnum = find(tracks2018.all.atype==-1); % ����������
CElt = tracks2018.all.lifetime(find(tracks2018.all.atype==1)); % ������������ lifetime
AElt =  tracks2018.all.lifetime(find(tracks2018.all.atype==-1)); % ��������������

% % ֱ��ͳ��ͼ(fig1)
% set(gcf,'Position',[100 50 700 250]);
% hold on;
% edges=[1:1:50]
% h1 = histogram(CElt,edges,'facecolor',[0 0 0]);
% h2 = histogram(AElt,edges,'facecolor',[1 1 1]);
% legend({'CE','AE' });
% box on
% xlabel('��������/��','fontweight','bold');
% ylabel('����','fontweight','bold');
% xlim([0,50])
% set(gca,'xtick',[0:5:50])
% % a=tabulate(tracks2018.all.lifetime(:)); % ͳ��Ƶ��
clear h1 h2 
%% ����SLA����
file='..\data\originaldata\sla_20180701_1231_lon_144.5_165.5_lat_37.5_50.5.nc';
ncdisp(file)
SLA2018.all.latitude = ncread(file,'latitude');
SLA2018.all.longitude = ncread(file,'longitude');
SLA2018.all.u = ncread(file,'ugos'); % ���� m/s
SLA2018.all.v = ncread(file,'vgos');
SLA2018.all.sla = ncread(file,'sla'); % ������߶��쳣ֵ
SLA2018.all.time = ncread(file,'time'); % days since 1950-01-01 00:00:00
SLA2018.all.day = SLA2018.all.time+datenum(1950,1,1) - datenum(2018,7,1)+1;% ƥ����ʦ������������
clear file

%% ƥ��SLA���ݺ�u v�������������ĵ�ͱ߽�
for i = 1:length(tracks2018.all.day)
    thiseddyday = cell2mat(tracks2018.all.day(i)); % ����
    thiseddylon = cell2mat(tracks2018.all.lon(i)); % ����
    thiseddylat = cell2mat(tracks2018.all.lat(i)); % γ��
    thiseddyedge = tracks2018.all.shapes(i); % �߽���Ϣ
    t1=[];
    t2=[];
    t3=[];
    t4=[];
    distance_ab=[];
    for j = 1:length(thiseddyday)
        d = thiseddyday(j);
        a1 = thiseddylon(j);
        a2 = thiseddylat(j);
        b1 = thiseddyedge{1}{j}(1,:); % ����
        b2 = thiseddyedge{1}{j}(2,:); % γ��
        distance_ab(j) = mean(distance(a1,a2,b1,b2)/180.*pi.*6371);% �����뾶�Ķ����������߽��ϸ��㵽���ĵľ���ƽ��ֵ��
        [lon2,lat2] = meshgrid(SLA2018.all.longitude,SLA2018.all.latitude);
        eu = interp2(lon2,lat2,SLA2018.all.u(:,:,d)',b1,b2,'spline'); % ƥ��ı߽�����
        ev = interp2(lon2,lat2,SLA2018.all.v(:,:,d)',b1,b2,'spline');
        csla =interp2(lon2,lat2,SLA2018.all.sla(:,:,d)',a1,a2,'spline'); % ƥ�������SLA
        esla =interp2(lon2,lat2,SLA2018.all.sla(:,:,d)',b1,b2,'spline'); % ƥ��ı߽�SLA
        t1(j) = csla;
        t2(j) = mean(esla); % �߽�sla��ֵ
        t3(j) = mean(esla)-csla; % ��� SLAƫ��
        % ��������ж�
        [lonm,latm]=meshgrid(b1,b2);
        [um,vm]=meshgrid(eu,ev);
        curlz=curlz_atmos(lonm,latm,um,vm);
        t4(j)=mean(curlz(:),'omitnan'); % �߽�����жȾ�ֵ
        
    end
    tracks2018.all.radius{i} = distance_ab; % �뾶
    tracks2018.all.centersla{i} = t1;
    tracks2018.all.edgesla{i} = t2;
    tracks2018.all.Resla{i} = t3; % ��� С��0��ʾ���������½�����ů������ɫ AE������0��ʾ���� CE ��ɫ
    tracks2018.all.Recurlz{i} = t4;  % ����ж�
    tracks2018.all.mappedday{i} = mapminmax(cell2mat(tracks2018.all.day(i))', 0, 1); % ��׼������������
end
clear t1 t2 t3 t4 a1 a2 AElt CElt b1 b2 d distance_ab 
clear thiseddyday thiseddyedge thiseddylat thiseddylon i j
clear um vm  eu ev flag lon2 lat2 latm lonm curlz esla csla

%% �������ʵ��ݱ� ���뾶�����������ж��ڱ�׼�����������ڵ�ƽ���ݱ������
%�������ɸ�����
ceer=[]; % eddy radius
ceec=[]; % eddy curl
cees=[]; % eddy relatively SLA
ceed=[]; % ��׼������������
aeer=[]; % �� radius
aeec=[]; % �� curl
aees=[]; % �� relatively SLA
aeed=[]; % �� ��׼������������
scday=[];saday=[];
sclat=[];sclon=[];salat=[];salon=[];
for i = 1:length(tracks2018.all.day)
    if tracks2018.all.atype(i)==1
        ceed =[ceed,cell2mat(tracks2018.all.mappedday(i))];% ��׼������������
        ceer =[ceer,cell2mat(tracks2018.all.radius(i))]; % �����뾶
        cees =[cees,cell2mat(tracks2018.all.Resla(i))]; % ����SLAƫ��
        ceec =[ceec,cell2mat(tracks2018.all.Recurlz(i))]; %��������ж�
        sclat=[sclat,cell2mat(tracks2018.all.lat(i))'];
        sclon=[sclon,cell2mat(tracks2018.all.lon(i))'];
        scday=[scday,cell2mat(tracks2018.all.day(i))'];
    elseif  tracks2018.all.atype(i)==-1
        aeed =[aeed,cell2mat(tracks2018.all.mappedday(i))];
        aeer =[aeer,cell2mat(tracks2018.all.radius(i))]; % �����뾶
        aees =[aees,cell2mat(tracks2018.all.Resla(i))]; % ����SLAƫ��
        aeec =[aeec,cell2mat(tracks2018.all.Recurlz(i))]; %��������ж�
        salat=[salat,cell2mat(tracks2018.all.lat(i))'];
        salon=[salon,cell2mat(tracks2018.all.lon(i))'];
        saday=[saday,cell2mat(tracks2018.all.day(i))'];
    end
end

CEpos=[ceed',ceer',cees',ceec',sclat',sclon']; % ��׼������������ �뾶 ��� ����ж� γ�� ����
AEpos=[aeed',aeer',aees',aeec',salat',salon'];

%  ����׼�����������ڵĴ�С����������ͬ���ڵ�������ȡƽ��ֵ
[bc mc nc]=unique(CEpos(:,1));
for ii=1:length(mc)
    bc(ii,2)=mean(CEpos(nc==ii,2)); % �뾶
    bc(ii,3)=mean(CEpos(nc==ii,3)); % ���
    bc(ii,4)=mean(CEpos(nc==ii,4)); % ����ж�
end

[ba ma na]=unique(AEpos(:,1));
for ii=1:length(ma)
    ba(ii,2)=mean(AEpos(na==ii,2));
    ba(ii,3)=mean(AEpos(na==ii,3));
    ba(ii,4)=mean(AEpos(na==ii,4));
end
mapedday2=linspace(0,1,20); % �������е�������[0, 1]������ƽ���ֳ�20��ʱ���
% ������������������ֵ����20��ʱ�����
cenewradius =abs(interp1(bc(:,1),bc(:,2),mapedday2));
aenewradius =abs(interp1(ba(:,1),ba(:,2),mapedday2));
cenewResla =abs(interp1(bc(:,1),bc(:,3),mapedday2));
aenewResla =abs(interp1(ba(:,1),ba(:,3),mapedday2));
cenewRecurlz =abs(interp1(bc(:,1),bc(:,4),mapedday2));
aenewRecurlz =abs(interp1(ba(:,1),ba(:,4),mapedday2));

clear aeec aeed aeer aees bc ba ceec ceed ceed ceer cees salat salon
clear i ii j mc ma nc na sclat sclon saday scday CEnum AEnum

% ��ͼ (fig2)
set(gcf,'Position',[100 20 1000 700]);
ylabletxt={'�뾶/km','\deltaSLA/m','����ж�'};
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
    xlabel('��׼������������','fontweight','bold');
    ylabel(ylabletxt{i},'fontweight','bold');
    xli=xlim;
    yli=ylim;
    text(xli(1)+0.03,yli(1)+0.9*(yli(2)-yli(1)),texttxt{i})
end
legend({'CE','AE' });
clear cenewResla aenewResla cenewRecurlz aenewRecurlz aenewradius cenewradius 
clear i k mapedday2 xli yli y1 y2 y3 ylabletxt texttxt

%% �������뾶������ĵ���ֲ�ͼ(fig3)
% set(gcf,'Position',[100 50 900 500]);
% set(gcf,'color','white');
% for i=1:6
%     % CEpos=[ceed',ceer',cees',ceec',sclat',sclon']; % ��׼������������ �뾶 ��� ����ж� γ�� ����
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
%         v=abs(d(:,2)); %�뾶
%     elseif i==3|i==4
%         v=abs(d(:,3)); %���
%     else
%         v=abs(d(:,4)); %����ж�
%     end
%     
%     [xq,yq] = meshgrid(140:1:170, 36:1:55); %���о��������ֳ�1���1������
%     vq = griddata(x,y,v,xq,yq); % , ������������������ƶ�������������ƽ���뾶�����������ж�
%     m_proj('equidist', 'lon',[146 164],'lat',[39 49], 'aspect', 0.5);
%     m_pcolor(xq,yq,vq);
%     m_coast('patch', [.7 .7 .7])
%     m_coast('linewidth', 1, 'color', 'b');%����������
%     m_grid('box','fancy');%��Ӹ���;
%     colorbar
%     colormap(jet);
%     hold on
%     xlabel('����','FontWeight','bold')
%     ylabel('γ��','FontWeight','bold')
%     % caxis([0 0.4])
%     h = colorbar;
% end
clear A ans C d h i ia ic v vq x xq y yq

%% ���뾶�����ܶ�ͼ (fig4)
% figure
% set(gcf,'Position',[100 50 700 250]);
% % CEpos=[ceed',ceer',cees',ceec',sclat',sclon']; % ��׼������������ �뾶 ��� ����ж� γ�� ����
% % AEpos=[aeed',aeer',aees',aeec',salat',salon'];
% %����
% y1 = CEpos(:,2);
% y1min=min(CEpos(:,2));
% y1max=max(CEpos(:,2));
% %������
% y2 = AEpos(:,2);
% y2min=min(AEpos(:,2));
% y2max=max(AEpos(:,2));
% x=linspace(y2min,y2max,20); % �������С����ֳ�20���ȷֵ�(19�ȷ�),
% yy1=hist(y1,x); % �����������ĸ���
% yy1=yy1/length(y1); % �����������ĸ���
% yy2=hist(y2,x);
% yy2=yy2/length(y2);
% % bar(x,yy);%���������ܶȷֲ�ͼ
% scatter(x,yy1,25,'b','LineWidth',1.5);
% hold on;
% scatter(x,yy2,25,'r','LineWidth',1.5);
% xlabel('�뾶/km','fontweight','bold');
% ylabel('�����ܶ�','fontweight','bold');
% set(gca,'ygrid','on');
% box on
% legend({'CE','AE' });
clear x y1 y1max y1min y2 y2max y2min yy1 yy2

%% ���������������ж���뾶�Ĺ�ϵ 
% ���հ뾶����
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
mapedrad2=linspace(10,110,20);%���뾶��С������5km���ֳ����ɸ�����,
cenewResla =abs(interp1(bc(:,1),bc(:,2),mapedrad2));
aenewResla =abs(interp1(ba(:,1),ba(:,2),mapedrad2));
cenewRecurlz =abs(interp1(bc(:,1),bc(:,3),mapedrad2));
aenewRecurlz =abs(interp1(ba(:,1),ba(:,3),mapedrad2));

% % ��ͼ(fig5)
% figure
% set(gcf,'Position',[100 50 700 400]);
% ylabletxt={'\deltaSLA/m','����ж�'};
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
% xlabel('�뾶/km','fontweight','bold');
% ylabel(ylabletxt{i},'fontweight','bold');
% xli=xlim;
% yli=ylim;
% text(xli(1)+0.2,yli(1)+0.9*(yli(2)-yli(1)),texttxt{i})
% end
% legend({'CE','AE' });
clear cenewResla aenewResla cenewRecurlz aenewRecurlz aenewradius cenewradius 
clear ba bc i ii k ma mc nc texttxt xli y1 y2 y3 ylabletxt yli na mapedrad2


%%  �о������ڸ��к�ƽ�涯���߶��쳣(fig6)
% j=0
% figure('position',[333,44,600,300]);
% for i=1:length(tracks2018.all.lifetime)
%     lt=tracks2018.all.lifetime(i)
%     if lt>=56&lt<63 % ��
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