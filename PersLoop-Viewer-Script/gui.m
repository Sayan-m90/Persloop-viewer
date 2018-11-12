filename = 'botijo_002_pers';
fid = fopen(filename,'r');
one = [];
count = 2;
color = [];
% cinf = [];
C=cell(1,3); C{1} = '0';
while(C{1}~='2')
    C = strsplit(fgetl(fid));
if (C{1}=='1')
    if(strcmp(C{3},'inf'))
%         oneinf=[oneinf;[str2num(C{2}),-10]];
%         cinf=[cinf;'black'];
%         color = c;
    else
        one=[one;[str2num(C{2}),str2num(C{3})]];
        c=[color;'b'];
        color = c;
    end
end
end
% hold on

[v1,f1] = read_off('botijo.off');
          

figure,h = barh(one,'stacked','EdgeColor',[1 1 1]);
fig=gcf;
set(h,{'FaceColor'},{'w';'b'});
btn = uicontrol('Style', 'pushbutton', 'String', 'Filtration',...
        'Position', [20 20 50 20],...
        'Callback', 'cla');
% c= get(0,'children');
% set(c,'DefaultFigurePosition', [0 0 0 0]);
% set(h, 'Position', [0 0]);
% movegui('east');    
%     fig=gcf;
% set(h, 'DefaultFigurePosition', [-10 -10 500 300]);    
while(true)
    display(btn);
    [xlist,ylist]=ginput(1);
    display(xlist)
    if xlist<=0
        continue;
    end
    lp = strcat(filename(1:end-5),'floops/');
    loopoff = strcat(strcat(lp,num2str(one(uint8(ylist),1))),'nc.off');
    display(loopoff);
    [v2,f2] = read_off(loopoff);
%     clf(2,'reset');
% set(1, 'DefaultFigurePosition', [100 100 500 300]);    
    plot_mesh2(v1, f1,[],v2,f2,[]);
%     display('out');
%     close(fig);
%     count=count+1;
    figure,h = barh(one,'stacked','EdgeColor',[1 1 1]);
%     movegui('east');
    fig=gcf;
    set(h,{'FaceColor'},{'w';'b'});
    % display(btn);
    % one(uint8(ylist):1)+;
end