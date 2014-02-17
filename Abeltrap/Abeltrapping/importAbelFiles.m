function AbelInfo = importAbelFiles(pathToData,filenameHeaderDat,filenameHeaderTXT,filenamePhotonDat,filenameFeedbackDat)
%function ABELinfo = importAbelFiles(pathToData,filenameHeaderDat,filenameHeaderTXT,filenamePhotonDat,filenameFeedbackDat)
%Import the parameters of an ABEL experiment and store them in the AbelinfoStruct
%The intention is to read multiple experiments simultaneously and make use
%of the cluster.
%
%input arguments:
%pathToData is the path to the stored files
%filenameHeaderDat, filenameHeaderTxt, filenamePhotonDat and
%filenamefeedbackDat are the experiment datafiles of one record
%
%return struct:
%ABELinfo is a struct with fields:
%pathToData, filenameHeaderDat, filenameHeaderTxt, filenamePhotonDat and,binspercycle,binlength,EODx,EODy,beamwaist,EODscale,diffconst,mobility,maxV,mobilityon
%Read the Parameter of the experiment   

%Get the handlers of the headerfiles
headerfile = fopen(fullfile(pathToData,filenameHeaderDat),'r');
textheaderfile = fopen(fullfile(pathToData,filenameHeaderTXT),'rt');
textheader=textscan(textheaderfile,'%s',49);
header = textheader{1};
tempbinlength =  header{9};
binlength=str2double(strrep(tempbinlength,',','.'));
tempbinspercycle = header{17};
binspercycle=str2double(strrep(tempbinspercycle,',','.'));
tempbeamwaist = header{28};
beamwaist=str2double(strrep(tempbeamwaist,',','.'));
tempEODscale = header{32};
EODscale=str2double(strrep(tempEODscale,',','.'));
tempdiffconst = header{36};
diffconst=str2double(strrep(tempdiffconst,',','.'));
tempmobility = header{39};
mobility=str2double(strrep(tempmobility,',','.'));
tempmaxV = header{43};
maxV=str2num(strrep(tempmaxV,',','.'));
mobilityon = header{48};
if(strcmp(mobilityon, 'FALSE'))
    mobilityon = false;
else
    mobilityon = true;
end
EODx = fread(headerfile,binspercycle,'double');
EODy = fread(headerfile,binspercycle,'double');
fclose(headerfile);
fclose(textheaderfile);
%Save the Imported Info
AbelInfo = struct('pathToData',pathToData,'headerDAT',filenameHeaderDat,'headerTXT',filenameHeaderTXT,'photonDat',filenamePhotonDat,'feedbackDat',filenameFeedbackDat,'binspercycle',binspercycle,'binlength',binlength,'EODx',EODx,'EODy',EODy,'beamwaist',beamwaist,'EODscale',EODscale,'diffconst',diffconst,'mobility',mobility,'maxV',maxV,'mobilityon',mobilityon);
