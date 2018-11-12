[filename,folder] = uigetfile({'*.off';'*.png';'*.jpg'},'Enter barcode file');
fid = fopen(strcat(strcat(folder,filename(1:end-4)),'_pers'));
ftx = fopen(strcat(strcat(folder,filename(1:end-4)),'f.txt'));
setenv('A',strcat(strcat(folder,filename(1:end-4)),'_pers'));
setenv('B',strcat(strcat(folder,filename(1:end-4)),'f.txt'));
!./persloop -s $A -f $B
offfile = fopen(strcat(folder,filename),'r');
picflag = 0;
if(strcmp(filename(end-3:end),'.off')==false)
    picflag = 1;
end
one = [];
count = 3;
color = [];
% cinf = [];
C=cell(1,3); C{1} = '0';
% if picflag == 0;
while(C{1}~='2')
%     display(fgetl(fid))
    C = strsplit(fgetl(fid));
if (C{1}=='1')
    if(~strcmp(C{3},'inf'))
        one=[one;[str2num(C{2}),str2num(C{3})]];
    end
end
end
% hold on
one(end+1,:)=NaN;
% end
% text(10,20,'Explanation');
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',size(one,1),'Value',size(one,1),...
        'Position', [50 20 150 20],'SliderStep',[.1 1]); 
toplot={'OFF','PNTS','FILT'};
pop = uicontrol('Tag','popupmenu','Style', 'popup','Position',[250 -10 100 50],...
'string',toplot);
% offfile = strcat(filename(1:end-5),'color.off');
% if (picflag == 0)
display(offfile);
if (picflag == 0)
[v1,f1,f2,f3,~] = read_off(strcat(folder,filename));

% display(size(f2));
% display(size(f3));
f = vertcat(f2,f3);
end
% f1=horzcat(f,f1);
% exit(0);
% load(strcat(filename(1:end-5),'var'));
% [vertex,fTri,fTri_time,fEdge,fEdge_time] = read_filtration(strcat(filename(1:end-5),'f.txt'));

ylist = -1;
% figure;
while(ishandle(gcf))
    set(gcf, 'Color', [1.0,1.0,1.0]);
    delete(subplot(1,2,2));
    subplot(1,2,2);    
    data = one(1:sld.Value-1,:);
    data(end+1,:)=NaN;
    b={barh(data([1 end],:),'stacked')};
    hold on
    for k=1:(size(data,1)-1)
        b{k}=barh([k-1 k],data([end k],:),'stacked','EdgeColor',[0.1,0.1,0.1]);
    end

    for k=1:(size(data,1)-1)
            set(b{k}(1),'FaceColor',[0.1,0.1,0.1])
            set(b{k}(2),'FaceColor',[0.8500, 0.3250, 0.0980])        
    end

    set(gcf, 'Color',[.1,.1,.1]);%[0.1,0.1,0.1]);%[1,1,1]);% 
      set(subplot(1,2,2), 'Color',[.1,.1,.1]);% [0.1,0.1,0.1]);%[1,1,1]);%,
      text(10000,2020,'Explanation','Color','r');
%     if ylist>=0
%         set(b{uint8(ylist)}(2),'FaceColor','r')  
%     end
    hold off;
%     text(-10,120,'Explanation','Color','b');
    [xlist,ylist]=myginput(1,'crosshair');
    display(uint8(ylist));
    lp = strcat(filename(1:end-4),'floops\');
    loopoff = strcat(strcat(lp,num2str(data(uint8(ylist),1))),'comb.off');
%     display(loopoff);
    delete(subplot(1,2,1));
    subplot(1,2,1);
%     text(-10,120,'Explanation','Color','w');
    hold on;
    display(loopoff);
    [v21,f21,~,~,c21] = read_off(strcat(folder,loopoff));%'mc_0018_tree_outputfloops/47137comb.off');%
%     [v22,f22,c22] = read_off('mc_0018_tree_outputfloops/48809comb.off');
%     [v23,f23,c23] = read_off('mc_0018_tree_outputfloops/46542comb.off');
    v2 = [v21];%,v22,v23];
    f2 = [f21];%,f22,f23];
    c2 = [c21];%,c22,c23];
    E = get(pop,{'value'});
    v22=[];
    f22=[];
    c22=[];
    v23=[];
    f23=[];
    c23=[];
    if picflag ==1
        imshow(strcat(folder,filename));
        hold on;
        XX = [];
        YY = [];
        v21 = v21';
        for i=1:size(f21,2)
            XX = [XX v21(f21(1,i),1) v21(f21(2,i),1) NaN];
            YY = [YY v21(f21(1,i),2) v21(f21(2,i),2) NaN];
%             ZZ = [ZZ vertex(face_edge(1,i),3) vertex(face_edge(2,i),3) NaN];
        end
        if(size(f21,1)~=0)
                line(YY,XX,'Color',[0.9500, 0.1250, 0.0980],'LineWidth',2);
        end
        hold off;

    elseif E{1,1}==1
%         plot_mesh2(v1, f1,[],v21,f21,c21,v22,f22,c22,v23,f23,c23);
        plot_mesh2(v1, f1,[],v2,f2,c2,v22,f22,c22,v23,f23,c23);
    elseif E{1,1}==2
        fz = f1(:,1);
%         plot_mesh2(v1, fz,[],v2,f2,c2);
    plot_mesh2(v1, fz,[],v21,f21,c21,v22,f22,c22,v23,f23,c23);
    else
%         [vertex,fTri,fEdge] = read_filtration('botijo_002f.txt',one(uint8(ylist),1));
        B = fEdge_time<data(uint8(ylist),1);
        k = find(B,1,'last');
        fE = fEdge(:,1:k);
        B = fTri_time < data(uint8(ylist),1);
        k = find(B,1,'last');
        fT=fTri(:,1:k);
        plot_mesh2(v1, fT, fE,v22,f22,c22,v21,f21,c21,v23,f23,c23);
%         plot_mesh2(v1, fT, fE, v2, f2, c2);
    end
    hold off;
    if ylist>=0
        set(b{uint8(ylist)}(2),'FaceColor',[0.4940, 0.1840, 0.5560])
    end
% text(-10,2000,offfile,'Color','w','Interpreter', 'none');
    pause;

%     count=get(gcf,'CurrenKey');
end