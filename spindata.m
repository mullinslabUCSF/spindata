function [nspeed,eqapp1,eqapp2,eqapp3,xm,xb] = spindata(datadir,outdir,fextens,fnstart)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R. Dyche Mullins 
% Updated: 11/10/2020
%
% This function starts the process of analyzing sedimentation equilibrium
% data. The input string 'datadir' points to a folder containing raw
% equilibrium sedimentation data, and 'fextens' contains the filename 
% extension of the data to be analyzed (e.g."RA1", "RA2", etc). A GUI 
% enables the user to choose two sets of boundary points for each of the 
% three channels in the centerpiece. The first set of boundary points is 
% used to parse the data for analysis and should be chosen conservatively. 
% For example, these boundary points should exclude noisy points near the 
% meniscus and any data points with absorbance >0.9 (the limit of the 
% linear range of the detector). The second set of boundary points is 
% used by other functions (e.g. multi_global_kb2()) to estimate the total 
% amount of protein in the cell and should reflect, as accurately as 
% possible, the actual edges of the channels. These boundary points 
% should always lie outside of the first set of boundary points. 
% 
% In addition, this function detects the total number of rotor speeds 
% used in the experiment; finds the final absorbance trace collected at 
% each speed; and plots the ‘distance from equilibrium’ for each data 
% set, measured as the mean-square difference from the final data set 
% collected at that speed. 
%
% Note: this function assumes that you are using 3-channel centerpieces and
% collecting data at three different speeds. If you are using a different
% number of speeds you will have to modify the code for displaying the
% plots (the problem comes around line 175 and is an easy fix)
%
% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% datadir - string with the name of the directory that contains all the
% outdir - string with the name of the output file directory
% fextens - string with file extension for the proper cell/modality
% fnstart - string with a filename tag for the output data files
% to be analyzed
% OUTPUT VALUES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nspeed -- number of speeds represented in the datasets
% eqapp1, eqapp2, eqapp3 --approaches to equilibrium
% xm -- vector with estimated meniscus positions
% xb -- vector with estimated cell base positions

[nlist,headers,numfiles] = get_datafile_names(datadir,fextens);
% nlist - list of filenames with the appropriate suffix
% headers - header information from the data files
% numfiles - number of data files

speed = zeros(numfiles,1);
temp = zeros(numfiles,1);

% set up a cell array to hold ALL of the data for quality control and
% analysis. positions 1 and 2 are the x and y data. positions 3 and 4 are
% the speed (in RPM) and temperature (in degrees C)
megillah = cell(numfiles,5);
% here is a string to hold the first line of the file header for use later
% when we write out the output files
strhead1 = '__';
% and here is a cell array to hold the specific header info for all files
strhead2 = cell(numfiles,1);

% pull out all the data and do a preliminary plot
figure(1)
hold on
for k=1:numfiles
    [x,y,r,s,t,strhead1,strhead2{k}]=getdata(nlist{k},datadir);
    speed(k) = s;
    temp(k) = t;
    % plot all the data sets on top of each other
    plot(x,y);
    % load up the megillah cell array
    megillah{k,1} = x;
    megillah{k,2} = y;
    megillah{k,3} = r;
    megillah{k,4} = s;
    megillah{k,5} = t;
end
% set the axes for the data plot
axis([5.6 7.4 -0.1 1.2]);

% get endpoints for parsing the three separate channels from the user
xf = zeros(3,1);    % position of first data point in parsed set
xl = zeros(3,1);    % position of last data point in parsed set
for k=1:3
    [xf(k),y] = ginput(1);
    plot(xf(k),y,'ro');
    [xl(k),y] = ginput(1);
    plot(xl(k),y,'ro');
end

% get endpoints for the three separate channels from the user
xm = zeros(3,1);    % user estimated position of meniscus
xb = zeros(3,1);    % user estimated position of the channel bottom
for k=1:3
    [xm(k),y] = ginput(1);
    plot(xm(k),y,'k^','MarkerFaceColor','k');
    [xb(k),y] = ginput(1);
    plot(xb(k),y,'k^','MarkerFaceColor','k');
end
hold off

% set up cell arrays to hold data from each of the three cells
% first column is the x data; second column is the y data; third column is
% the std deviation
subcell1 = cell(numfiles,3);
subcell2 = cell(numfiles,3);
subcell3 = cell(numfiles,3);
% go parse the full data sets into the three sample cells
for k=1:numfiles
    [x1,y1,r1,x2,y2,r2,x3,y3,r3] = subreddit(megillah{k,1},megillah{k,2},megillah{k,3},xf,xl);
    % subreddit uses xm and xb to window the data
    % load up the subcells
    subcell1{k,1}=x1;
    subcell1{k,2}=y1;
    subcell1{k,3}=r1;
    subcell2{k,1}=x2;
    subcell2{k,2}=y2;
    subcell2{k,3}=r2;
    subcell3{k,1}=x3;
    subcell3{k,2}=y3;
    subcell3{k,3}=r3;
end

% replot the truncated data
figure(2)
subplot(2,1,1)
hold on
for k=1:numfiles
    plot(subcell1{k,1},subcell1{k,2});
    plot(subcell2{k,1},subcell2{k,2});
    plot(subcell3{k,1},subcell3{k,2});
end
title('all parsed curves');
    xlabel('radius (cm)');
    ylabel('absorbance (A.U.)');
hold off

[nspeed,splist,endpts] = eqfind(speed);
% nspeed - number of different speeds
% splist - list of indices for each of the centrifuge speeds
% endpts - indices of the final runs at each speed 

%figure(3)
subplot(2,1,2)
hold on
for k=1:nspeed
    plot(subcell1{endpts(k),1},subcell1{endpts(k),2});
    plot(subcell2{endpts(k),1},subcell2{endpts(k),2});
    plot(subcell3{endpts(k),1},subcell3{endpts(k),2});
end
title('parsed data at equilibrium');
    xlabel('radius (cm)');
    ylabel('absorbance (A.U.)');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the approach to equilibrium for all of the speeds
% equapp1-3 are cell arrays to hold the approaches to equilibrium
eqapp1 = cell(nspeed,2);
eqapp2 = cell(nspeed,2);
eqapp3 = cell(nspeed,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spintest: (1) regularizes x-axes and recalculates data by interpolation;
% (2)finds all centrifuge speeds used in the experiment; and (3) computes
% the MSD between the last curve at each speed and all preceding curves
eqapp1 = spintest(numfiles,subcell1,nspeed,endpts);
eqapp2 = spintest(numfiles,subcell2,nspeed,endpts);
eqapp3 = spintest(numfiles,subcell3,nspeed,endpts);

% plot the approach to equilibrium data for all speeds and channels on one
% composite figure, together with teh speeds
figure(4)
hold on
subplot(2,2,1);
plot([1:numfiles],speed);
axis([0 numfiles 0 30000]);
title('speed sequence');
    xlabel('scan');
    ylabel('speed (RPM)');
    hold off
for k=1:nspeed
    subplot(2,2,k+1);   % NOTE that this line will screw up if nspeed~=3
    hold on
    plot(eqapp1{k,1},eqapp1{k,2});
    plot(eqapp1{k,1},eqapp1{k,2},'ok');
    plot(eqapp2{k,1},eqapp2{k,2});
    plot(eqapp2{k,1},eqapp2{k,2},'ok');
    plot(eqapp3{k,1},eqapp3{k,2});
    plot(eqapp3{k,1},eqapp3{k,2},'ok');
    sttl = [num2str(speed(endpts(k))) ' RPM data'];
    title(sttl);
    xlabel('scan');
    ylabel('mean square difference (AU)');
    hold off
end
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finally, save the parsed data into text files that can be read by
% Winnonln or other analysis software
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the return directory to the current one
ret_dir = pwd;
% go to the output directory
cd(outdir)
for k=1:nspeed
    % convert speed to text (spstr)
    spstr = num2str(speed(endpts(k)));
    % format for headers
    hdrfmt = '%s \n';
    % format for saving the data
    datfmt = '   %7.5f    %7.5f    %7.5f \n';
    % assemble good names for all the files
    ch1name = [fnstart spstr '_CH1.txt'];
    ch2name = [fnstart spstr '_CH2.txt'];
    ch3name = [fnstart spstr '_CH3.txt'];
    
    % open and write the first data file
    fileID = fopen(ch1name,'w');
    % write the two header lines
    fprintf(fileID,hdrfmt,strhead1);
    fprintf(fileID,hdrfmt,strhead2{endpts(k)});
    % interleave the data so that it will work with fprintf
    a=subcell1{endpts(k),1};
    b=subcell1{endpts(k),2};
    c=subcell1{endpts(k),3};
    row_interleave = reshape([a(:) b(:) c(:)]',3*size(a,1), []);
    % write data to the file
    fprintf(fileID,datfmt,row_interleave);
    fclose(fileID);
    
    % open and write the second data file
    fileID = fopen(ch2name,'w');
    % write the two header lines
    fprintf(fileID,hdrfmt,strhead1);
    fprintf(fileID,hdrfmt,strhead2{endpts(k)});    
    % interleave the data so that it will work with fprintf
    a=subcell2{endpts(k),1};
    b=subcell2{endpts(k),2};
    c=subcell2{endpts(k),3};
    row_interleave = reshape([a(:) b(:) c(:)]',3*size(a,1), []);
    % write data to the file
    fprintf(fileID,datfmt,row_interleave);
    fclose(fileID);
    
    % open and write the third data file
    fileID = fopen(ch3name,'w');
    % write the two header lines
    fprintf(fileID,hdrfmt,strhead1);
    fprintf(fileID,hdrfmt,strhead2{endpts(k)});    
    % interleave the data so that it will work with fprintf
    a=subcell3{endpts(k),1};
    b=subcell3{endpts(k),2};
    c=subcell3{endpts(k),3};
    row_interleave = reshape([a(:) b(:) c(:)]',3*size(a,1), []);
    % write data to the file
    fprintf(fileID,datfmt,row_interleave);
    fclose(fileID);
end
% go back to the original directory
cd(ret_dir)

return
end

function [namelist,fcontent,numfiles] = get_datafile_names(data_dir,extens)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_datafile_names - goes to a directory; identifies all the .rax files
% and loads their names into a cell string array. 
% data_dir - where to look for the data files
% extens - string containing the extension for the files of interest (e.g.
% 'RA1', 'RA2', etc.
% RDM - 10/23/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the return directory to the current one
ret_dir = pwd;
% go to the specified data directory
cd(data_dir)
% do a linux 'ls' command and store the results in 'txtfiles'. This
% makes txtfiles a list of strings corresponding to the file names.
extens = ['.' extens];
searstr = ['*' extens];    % first make a search string
raxfiles = ls(searstr);

% make sure raxfiles is a single row vector
raxfiles=raxfiles';
raxfiles=reshape(raxfiles,1,[]);

% from the list of all files in the directory identify the .rax files and 
% store their positions in 'ind'. 
ind = findstr(extens,raxfiles);     % find the .RAX filename positions
% count the number of files
numfiles = length(ind);

% create an indexable list of namess that can be used to read the data
% files or rearranged to construct new names, linked to he originals. 
% This list of name parts gets passed on to other functions
namelist = cell(length(ind),1);
% create a string containing all the filenames 
for i = 1:length(ind)
    % Convert each element of txtfiles to a simpler name that is just the 
    % last three digits of the original name (usually three numbers) plus 
    % the '.RAX' string. 
    namelist{i} = raxfiles(1,ind(i)-5:ind(i)+3); % works for XLI files
end
% in case the files are out of order, sort them by filename
namelist = sort(namelist);

% get info about the contents of the data files
% column 1 is the file names; column 2 is the list of speeds
fcontent = cell(length(ind),2);
for k=1:numfiles
    % open the data file
    fileID = fopen(namelist{k},'r');
    dum = fgetl(fileID);
    fcontent{k,1} = fgetl(fileID);
    fclose(fileID);
    
    % find the speed information
    % find all the spaces in the string
    SPind = findstr(char(32),fcontent{k,1}); % ascii code for space is 32
    dumspeed = fcontent{k,1}(SPind(3)+1:SPind(4)-1);
    fcontent{k,2} = str2num(dumspeed);
    
end

cd(ret_dir)

return
end

function [x,y,err,speed,temp,sh1,sh2] = getdata(datafile,datadir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open an analytical ultracentrifuge data file and copy out all the data
%
% x - x positional data
% y - absorbance at each position in the cell
% err - calculated standard deviation at each point
% speed - centrifuge speed (determined from header information
% sh1 - first line of the header
% sh2 - second line of the file header
%
% datafile - is the filename to open
% datadir - is the directory that holds all the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the return directory to the current one
ret_dir = pwd;
% go to the specified data directory
cd(datadir)

% get info about the contents of the data files
% column 1 is the file names; column 2 is the list of speeds

% open the data file
fileID = fopen(datafile,'r');
sh1 = fgetl(fileID);    % read the first line of the header
sh2 = fgetl(fileID);    % read the second line of the header
    
% find the speed information
% find all the spaces in the string
SPind = findstr(char(32),sh2); % ascii code for space is 32
dumparam = sh2(SPind(3)+1:SPind(4)-1);
speed = str2num(dumparam);
% the code above fails when the speed is below 10000 RPM so let's fix it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(speed)
    dumparam = sh2(SPind(4)+1:SPind(5)-1);
    speed = str2num(dumparam);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dumparam = sh2(SPind(2)+1:SPind(3)-1);
temp = str2num(dumparam);

% count the number of data points
fcontent = fgetl(fileID);
k = 0;
while(fcontent ~= -1)
    k=k+1;
    fcontent = fgetl(fileID);
end

% set up data arrays
numpts = k;
x = zeros(numpts,1);
y = zeros(numpts,1);
err = zeros(numpts,1);

fclose(fileID);
fileID = fopen(datafile,'r');
% skip header
dum = fgetl(fileID);
dum = fgetl(fileID);
for k=1:numpts
    fcontent = fgetl(fileID);
    A = sscanf(fcontent,'  %f  %f   %f');
    x(k) = A(1);
    y(k) = A(2);
    err(k) = A(3);
end

% figure()
% plot(x,y);
% cd(ret_dir)

cd(ret_dir)

end

function [x1,y1,r1,x2,y2,r2,x3,y3,r3] = subreddit(x,y,r,xm,xb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% subreddit uses the meniscus and bottom positions (xm and xb) supplied by
% the user in the parent algorithm to window the data in arrays x and y,
% along with the std dev in the array r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numpts = length(x);

% count the number of points in each window. the len_win array holds the
% lengths of the windows
win_len = zeros(3,1);
for k=1:numpts
    for j=1:3
        if x(k)>xm(j) && x(k)<xb(j)
            win_len(j)=win_len(j)+1;
        end
    end
end

% initialize the output vectors
x1 = zeros(win_len(1),1);
y1 = zeros(win_len(1),1);
r1 = zeros(win_len(1),1);
x2 = zeros(win_len(2),1);
y2 = zeros(win_len(2),1);
r2 = zeros(win_len(2),1);
x3 = zeros(win_len(3),1);
y3 = zeros(win_len(3),1);
r3 = zeros(win_len(3),1);
% now do the assigning
l1=0;   % loop indices
l2=0;
l3=0;
for k=1:numpts
    if x(k)>xm(1) && x(k)<xb(1)
        l1=l1+1;
        x1(l1) = x(k);
        y1(l1) = y(k);
        r1(l1) = r(k);
    elseif x(k)>xm(2) && x(k)<xb(2)
        l2=l2+1;
        x2(l2) = x(k);
        y2(l2) = y(k);
        r2(l2) = r(k);
    elseif x(k)>xm(3) && x(k)<xb(3)
        l3=l3+1;
        x3(l3) = x(k);
        y3(l3) = y(k);
        r3(l3) = r(k);
    end
end

return
end

function [numsp,splist,ends] = eqfind(speeds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the number of unique speeds in a run and identify the last points
% speeds - list of speeds in the experiment
%
% numsp - number of different speeds
% splist - list of the unique speeds used in the experiment
% ends - indices of the final runs at each speed 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

splist = unique(speeds);
numsp = length(splist);
% initialize splist
ends=zeros(numsp,1);
% populate splist
for k=1:numsp
    ends(k) = find(speeds==splist(k),1,'last');
end

return
end

function eqout= spintest(nscans,channel,numsp,ends)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spintest - computes the mean square error between the last speed and each
% of the preceding speeds
%
% nscans - total number of centrifuge scans
% channel - data for all the scans in one (edited) channel
% numsp - number of speeds in this experiment
% ends - vector containing the numbers of final runs at each speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect stats on the scan positions for all the data sets
% npts is the maximum number of data points in all the scans
% loball is the lowest x value in the scan range
% hiball is the highest x value in the scan range
npts = 0;
loball = channel{1,1}(1);
hiball = channel{1,1}(length(channel{1,1}));
for k=1:nscans
    lomid = min(min(channel{k,1}));
    loball = max([lomid; loball]);
    himid = max(max(channel{k,1}));
    hiball = min([himid; hiball]);
    npts = max(max([length(channel{k,1}); npts]));
end

% equalize the sampling and interpolate the data
% start by creating a standard x vector
dstd = (hiball - loball)/(npts-1);
xstd = [0:(npts-1)]';
xstd = dstd*xstd;
xstd = xstd + loball;

% copy the channel data into another cell array for interpolation
interpolate = cell(nscans,3);
% re-sample and interpolate the data
for k=1:nscans
    interpolate{k,1}=xstd;
    for j=1:npts
        if find(channel{k,1}==xstd(j))
            % what to do if you get lucky
            luckind = find(channel{k,1}==xstd(j));
            interpolate{k,2}(j)=channel{k,2}(luckind);
        else
            % we want to find the two points immediately adjacent to the
            % current data point. To do this we first split the data set
            % into two parts, one above and one below the current point
            lowind=find(channel{k,1}<=xstd(j));
            lowerX=channel{k,1}(lowind);
            lowerY=channel{k,2}(lowind);
            
            % now let's find the data beyond the current point
            higind=find(channel{k,1}>=xstd(j));
            upperX=channel{k,1}(higind);
            upperY=channel{k,2}(higind);
            
            % here are the adjacent x values
            x0 = lowerX(length(lowerX));
            x1 = upperX(1);
            % and here are the corresponding y values
            y0 = lowerY(length(lowerY));
            y1 = upperY(1);
            % calculate the fractional position of the data point
            fup = (x1-xstd(j))/(x1-x0);
            flo = 1 - fup;
            % calculate the y value by linear interpolation
            interpolate{k,2}(j) = fup*y0 + flo*y1;
        end
    end
end
% output cell array
eqout = cell(numsp,2);

% now compute the approach to equilibrium for each run
% figure();
% hold on
firstrun = 0;
for k=1:numsp
    % count the number of runs at this speed
    nrun = ends(k)-firstrun;
    % initialize the approach vectors
    xapp=[1:nrun-1]';
    yapp=zeros(nrun-1,1);
    for j=1:nrun-1
        % calculate the mean square error
        % this is where the regularization of the x-axis comes in handy
        delvect=interpolate{ends(k),2}-interpolate{firstrun+j,2};
        yapp(j)=delvect*delvect';
        yapp(j)=yapp(j)/npts;
    end
    firstrun = ends(k);
    eqout{k,1}=xapp;
    eqout{k,2}=yapp;
%     plot(xapp,yapp);
%     plot(xapp,yapp,'ok');
end
% hold off
return
end


