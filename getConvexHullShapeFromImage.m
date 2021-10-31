%% This function performs an analysis of the physical properties of a given image and provides: 
    %% 1. The number of vertices on the convex hull (i.e., CH-Shape)
    %% 2. The convex hull area (i.e., CH-Area)
    %% 3. The number of items (i.e., numerosity)
    %% 4. A table summarizing these results. 
%% Note 1: some of the sub-functions require Image Processing Toolbox, but there are other alternatives for it.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Note 2: don't forget to set a working directory in line 42
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Outline: 
%{
 % Load dir
 % Count the number of total images and create a new cell that will be store all of the data (pre-allocation). 
 % Loop over the directory 
    % Load image file 
    % Convert to a binary matrix
    % Optional: Divide the image according to the bisecting line
    % Test the number of objects (numerosity)
    % Check the num' of vertices on the convex hull 
    % Store: image name, left convex hull shape (CHvL), right convex hull shape (CHvR), left convex hull area (CHareaL), left convex hull area(CHareaL)
 % Save csv file - writetable(data, filename)
%}

function getConvexHullShapeFromImage

%% Dock the figure to the MATLAB window
set(0,'DefaultFigureWindowStyle','docked')
disp('DefaultFigureWindowStyle')

flagPlot = 1; %% put this on zero - no plot
flagBiseceted = 1;  %% put this on zero - avoid removing the line
flagClose = 1; 
if flagClose == 1
    close all; 
end

%% Set directory and image list - don't forget to change to the relevant directory 
% Don't forget to check the folder
disp('change to correct folder') % just a warning message
directory = ''; 
imgs = collectImgsFromDir(directory);

%% Preallocate cell for data
data = preAllocateCellStruct(length(imgs)); 

%% Loop over the folder 
for j=length(imgs):-1:1 % Efficient preallocation: run from the end of the folder backward to its strart/top
    [thisName, thisFile] = currImg(imgs, j, directory);  % collect the current image attributes 
    
    %% After reading the content of the folder run and check for errors, convert to a Binary image and display it
    try 
    img = imread(thisFile, 'jpg');  % try to read image
    mat = convert2Bin(img);  %'im2bw' requires Image Processing Toolbox.
    figure('visible','off'), imshow(img); 
    title(thisName) ; 
    catch
    end
    
    %% Divide the image into two halves along the bisecting line and assign and show the divided images
    if flagBiseceted
    
       [imgL, imgR, lSide, rSide, numLSide, numRSide, midrow, midcol ] = divImg(mat);
        
    end
    
    %% Default option - Show full image
    mainImg =  mat; 
    currentImg = figure('visible','on'), imshow(mainImg);
    numCurr = get(gcf,'Number');
    
    %% Alternative Option
    % imgL = imcrop(mat,[0 0  midcol size(mat,1) ]);
    % lSide = figure('visible','on'), imshow(imgL); 
    % numLSide = get(gcf,'Number');

    %% Count the number of objects in each of the halves
    [nLeft, nRight, sL, sR] = nCount(imgL, imgR, thisName); 
    dImgL = sL;  dImgR = sR;
    
    %% count the number of vertices on the convex hull (CH-Shape)
    [nVertL, nVertR] = verNumCH(dImgL, dImgR);
    
    %% Calculate CH-area 
    [dImgLx, dImgLy, dImgRx, dImgRy, kL, kR , vL, vR, kLnew] = calcCHareas(dImgL, dImgR);
    
    %% Plot CH's vertices centroids and fill them
    plotCHareas(flagPlot, lSide, dImgLx, dImgLy, kL , rSide, dImgRx, dImgRy, kR, kLnew, numLSide, numRSide, j, imgs );

    data(j+1,1:7) = {thisName, nLeft , nRight , nVertL , nVertR  vL, vR }  ; 
end

%% Convert the data-cell into a table. The first row is reserved for variable names
data = cell2table(data);

%% Write the table to a CSV file
ratio = input('Type in the number ratio   ');
filename = ['Testing ' ratio  '.csv'];
writetable(data,filename);
end

function imgs = collectImgsFromDir(directory);
%% This function sets the relevant dir and collect the images from there as a struct 
imgs = dir(directory); % Read the content of the dir and assign to a struct
imgs = imgs(3:end); % Remove the 1st and 2nd files. They're windows-junk named : "." and ".."
end

%% Subfunctions:

function data = preAllocateCellStruct( rowNum )  
%%  This function prealloctes a cell for later. 
%{
   Preallocate cell/struct for data:
    nROw= number of imgs+ 1 (for the colNames)    ,    nCol = {name, CHvL, CHvR, CHareaL, CHareaL}
%}
colNames = {'Name', 'nObkjectLeft', 'nObkjectRIght', 'nVerticesLeft', 'nVerticesRight', 'CHareaLeft', 'CHareaRight' } ; 
data = cell(rowNum+1, 7);
data(1,:) = colNames;  
end

function [thisName, thisFile]= currImg(imgs, j, directory)
%% This function collects the current img attributes 
%     close all; 
    thisName = imgs(j).name; % Read the img name
    thisFile = fullfile(directory, thisName); % Read the file name including the path
end

function mat = convert2Bin(img) 
%%  This function loads and converts the image data as zeros & ones reflecting black and white pixels
    mat=rgb2gray(img); % convert to grayscale
    mat= im2bw(mat); 
end

function [imgL, imgR, lSide, rSide, numLSide, numRSide, midrow, midcol ] = divImg(mat)
%%  This function divides the image into its two halves along the bisecting line according to the location of the mid pixel (in x-y cartesian coordinates) [midrow, midcol, mat] = treatMid(mat);
imgL = imcrop(mat,[0 0  midcol-1 size(mat,1) ]);
imgR = imcrop(mat,[midcol+2 0  size(mat,2) size(mat,1)  ]);

lSide = figure('visible','on'), imshow(imgL); 
numLSide = get(gcf,'Number');

rSide = figure('visible','on'), imshow(imgR);
numRSide = get(gcf,'Number');
end 

function [midrow, midcol, mat] = treatMid(mat)
%% This function locates and removes the bisecting line
%% first find the indexes of the middle
midrow = ceil(size(mat,1)/2);
midcol = ceil(size(mat,2)/2);
%% Now delete it so U will get rid of the bisecting line 
mat = deleteMid(mat);
end 

function mat = deleteMid(mat)
%% A function that deletes the mid-bisecting line 
[row,colOfLine] = find(all(mat == 1,1));
mat = [ mat(:,1:(colOfLine(1)-1)), mat(:,colOfLine(end)+1:size(mat,2))   ];
end

function [nLeft, nRight, sL, sR] = nCount(imgL, imgR, ~) 
%{
  This function count the number of objects in each of the halves
    p.s. Later add a subfunction that will check if this matches the name
    of the file that also contain the data about the objects this should be
    done for Reliability considirations
    % function nCompare
%}  
%% Left side
    %% fill the gaps
    imgL = imfill(imgL, 'holes');
    bw = bwconncomp(imgL,8); 
    nLeft  = bw.NumObjects;
    %% find center of vertices coordinates 
    sL= regionprops(bw,'centroid');
    sL =cell2mat(struct2cell(sL));
    sL = [sL(1:2:end-1)' , sL(2:2:end)'];
%% Right Side
    %%    fill the gaps 
    imgR = imfill(imgR, 'holes');
    cc = bwconncomp(imgR,8);
    nRight  = cc.NumObjects; 
    %% find center verices coordinates
    sR= regionprops(cc,'centroid');
    sR =cell2mat(struct2cell(sR));
    sR = [sR(1:2:end-1)' , sR(2:2:end)'];
%{
   For additional reading
        Recommended: https://www.youtube.com/watch?v=Gq7mp3G94ao
        https://www.mathworks.com/matlabcentral/answers/110287-how-to-count-the-number-of-object-present-in-binary-image
        https://www.mathworks.com/help/images/ref/bwconncomp.html
        https://blogs.mathworks.com/videos/2009/10/09/finding-the-area-inside-a-convex-hull/
%}
end 

function [nVertL, nVertR] = verNumCH(dImgL, dImgR)
%% This function counts the number of vertices on the CH 
VertL = convhulln(dImgL);
VertR = convhulln(dImgR);

nVertL = size(VertL,1); 
nVertR = size(VertR,1); 

 
end 


function  [dImgLx, dImgLy, dImgRx, dImgRy, kL, kR , vL, vR, kLnew] = calcCHareas(dImgL, dImgR)
%% This function calculates the CH areas
dImgLx = dImgL(:,1); dImgLy = dImgL(:,2); 
dImgRx = dImgR(:,1); dImgRy = dImgR(:,2); 
[kL , vL] = convhull(dImgLx,dImgLy);
[kR, vR] = convhull(dImgRx,dImgRy);

drawDotInd = randperm(size(kL(2:end),1));
kLnew = kL(2:end); 
kLnew(drawDotInd(1),:) = [];
end 

function plotCHareas(flagPlot, lSide, dImgLx, dImgLy, kL , rSide, dImgRx, dImgRy, kR, kLnew, numLSide, numRSide, j, imgs )
%% This function Plot CH's vertices centroids and fill them
if flagPlot ==1
    
    %% Plot the Left side 
    % figure('visible','on') = (lSide); 
    set(groot,'CurrentFigure',numLSide);
    hold on; 
    plot(dImgLx(kL),dImgLy(kL),'.');
    axis equal;
    fill ( dImgLx(kL),dImgLy(kL), 'r','facealpha', 0.3 ); 
    title(['Left half- ' imgs(j).name]);
    hold off;

    %% Plot the Right side 
    set(groot,'CurrentFigure',numRSide);
    hold on; 
    plot(dImgRx(kR),dImgRy(kR),'.');
    axis equal;
    fill( dImgRx(kR),dImgRy(kR), 'g','facealpha', 0.3 ); 
    title(['Right half- ' imgs(j).name]);
    hold off;
    
end

end 



