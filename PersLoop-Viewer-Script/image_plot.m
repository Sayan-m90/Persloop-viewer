function image_plot(I,Icomb)
% I = imread('input/train-vol2-2.png');
imshow(I);
hold on;
[vertex,face_edge,face2,face3,color] = read_off('input/527151comb.off');
XX = [];
YY = [];
vertex = vertex';
for i=1:size(face_edge,2)
            XX = [XX vertex(face_edge(1,i),1) vertex(face_edge(2,i),1) NaN];
            YY = [YY vertex(face_edge(1,i),2) vertex(face_edge(2,i),2) NaN];
%             ZZ = [ZZ vertex(face_edge(1,i),3) vertex(face_edge(2,i),3) NaN];
end
if(size(face_edge,1)~=0)
            line(YY,XX,'Color',[0.9500, 0.1250, 0.0980],'LineWidth',2);
end
end