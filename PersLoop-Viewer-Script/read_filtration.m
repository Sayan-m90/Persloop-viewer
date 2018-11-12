function [vertex,fTri,fTri_time,fEdge,fEdge_time] = read_filtration(filename)
% read_filtration - read data from filtration file.
%
%   [vertex,face] = read_filtration(filename);
%
%   'vTri' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'fTri' is a 'nb.face x 3' array specifying the connectivity of the triangles.
%   'fEdge' is a 'nb.face x 3' array specifying the connectivity of the edges.
%
fid = fopen(filename,'r');
fEdge = [];
fTri = [];
fEdge_time = [];
fTri_time = [];
if( fid==-1 )
    error('Can''t open the file.');
    return;
end
str = fgets(fid);
[a,str] = strtok(str); ndim = str2num(a);
[a,str] = strtok(str); nvert = str2num(a);
% display(nvert);
[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;

str=fgets(fid);
% display(str);
str=fgets(fid); [~,str] = strtok(str); timehere = str2num(str);
% display(str);
% input('num');
while (~feof(fid))
% display(timehere);

str=fgets(fid);
% display(str);
% input('num');
if(size(str,2)==0 || size(str,1)==0 || strcmp(str,'-1'))
    break
end
C=strsplit(str);%display(C);
% return;
% if(C{1}=='c')
%     continue
% end
if(numel(C)==4)
    a = str2num(C{2});
    b = str2num(C{3});
    fEdge = horzcat(fEdge,[a+1;b+1]);
    fEdge_time = horzcat(fEdge_time,timehere);
%     display(fEdge);
%     display(fEdge_time);
%     display(timehere);
%     input('hey');
elseif(numel(C)==5)
    a = str2num(C{2});
    b = str2num(C{3});
    c = str2num(C{4});
    fTri = horzcat(fTri,[a+1;b+1;c+1]);
    fTri_time = horzcat(fTri_time,timehere);
end
% return;
if(feof(fid))
    break;
end
str=fgets(fid);
% display(str);
% input('num');

if(size(str,2)==0 || size(str,1)==0 ||strcmp(str,'-1'))
    break
end
[~,str] = strtok(str); 
timehere = str2num(str);
end

end