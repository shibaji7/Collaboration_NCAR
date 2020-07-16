%This function makes a simple avi movie from a sequence of frames
%The user can control the display time of each frame. The movie is created
%in the same folder where this function is run.
%
%Usage:
%      Inputs:          
%             name: root name of the framse 
%         filetype: '.bmp'   ,  '.jpg'  , '.jpeg' , '.tif'  ,,,,
% number_of_frames: number of the frames
% display_time_of_frame: it is the display time of each frame  
%           
%           Output:
%  it creates a file named as Movie.avi in the same folder of this function
%           
%            Usage: Try the simple provided example
%   Movie_from_frames('image','.bmp',4,10)
%This function is written by :
%                             Nassim Khaled
%                             American University of Beirut
%                             Research Assistant and Masters Student
%Modified from the original


function Movie(display_time_frame)
mov = VideoWriter('raytrace.avi');
open(mov);
display_time_of_frame = display_time_frame;
myFolder ='.\raytrace-north';
filePattern = fullfile(myFolder, '*.png');
theFiles = dir(filePattern);
% for k= 1:length(theFiles)
% theFiles(k).img_no = str2num(theFiles(k).name((strfind(theFiles(k).name, '_')+1):(strfind(theFiles(k).name, '.')-1)));
% end
% cell_theFiles = struct2cell(theFiles);
% theFiles = sort({theFiles.img_no});
write_movie(mov, myFolder, theFiles,display_time_of_frame);
myFolder ='.\raytrace-west';
filePattern = fullfile(myFolder, '*.png');
theFiles = dir(filePattern);
write_movie(mov, myFolder, theFiles,display_time_of_frame);
myFolder ='.\raytrace-east';
filePattern = fullfile(myFolder, '*.png');
theFiles = dir(filePattern);
write_movie(mov, myFolder,theFiles,display_time_of_frame);
close all
close(mov);
end

function write_movie(mov, myFolder, theFiles,display_time_of_frame)
count=0;
for i=1 : length(theFiles)
    %name1=strcat(name,num2str(i),filetype);
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    a=imread(fullFileName);
    while count<display_time_of_frame
        count=count+1;
        imshow(a);
        F=getframe(gca);
        writeVideo(mov,F);
    end
    count=0;
end
end
