myDir = uigetdir;
myFiles = dir(fullfile(myDir,'*.png'));
% se = strel('disk',3);
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  display(baseFileName);
  J = imread(strcat('all-images/',baseFileName));
%   JK=imclose(J,se);
%   JK = bwmorph(JK,'skel',Inf);
%   imwrite(JK,strcat('all-images/',strcat(baseFileName(1:end-4),'m.png')));
%   delete(strcat('all-images/',baseFileName));
[h,w]=size(J);
numberOfTruePixels = sum(J(:));
fileID = fopen(strcat('all-images/',strcat(baseFileName(1:end-4),'.off')),'w');
fprintf(fileID,'OFF\n%d 0 0\n',numberOfTruePixels);
for i=1:h
    for j=1:w
        if(J(i,j)==1)
            fprintf(fileID,'%d %d 0\n',i,j);
        end
    end
end
fclose(fileID);
end
