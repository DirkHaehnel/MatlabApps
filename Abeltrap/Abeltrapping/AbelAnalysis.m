

clear all;
close all;
%Open the datafiles in a dialog thus nothing can go the wrong way
[filenameHeaderDat, pathToData] = uigetfile('*.dat', 'Pick an AbelTrap-File Header-Log-File in DAT format','\\jesrv\WG-Data\*');
if isequal(filenameHeaderDat,0)
    return
else
[filenameHeaderTXT, pathToData] = uigetfile('*.dat', 'Pick an AbelTrap-File Header-Log-File in TXT format','\\jesrv\WG-Data\*');
if isequal(filenameHeaderTXT,0)
    return
else
[filenamePhotonDat, pathToData] = uigetfile('*.dat', 'Pick an AbelTrap-File Photon-Log-File in DAT format','\\jesrv\WG-Data\*');
if isequal(filenamePhotonDat,0)
    return
else
[filenameFeedbackDat, pathToData] = uigetfile('*.dat', 'Pick an AbelTrap-File Feedback-Log-File in DAT format','\\jesrv\WG-Data\*');
if isequal(filenameFeedbackDat,0)
    return
else
%call datafiler
AbelInfo = importAbelFiles(pathToData,filenameHeaderDat,filenameHeaderTXT,filenamePhotonDat,filenameFeedbackDat);
%find trapped molecules
AbelData = traceAbelDataFiles(AbelInfo);

AbelTraping = findTrappedMolecules(AbelData.photonbins, AbelInfo.binspercycle, AbelInfo.binlength);





%Plot the whole thing
coarsegraintime = 10e-4; %seconds
binspercoarsegrain = AbelInfo.binspercycle*round(coarsegraintime/(AbelInfo.binspercycle*AbelInfo.binlength*1e-9));
    %Ensure an integer number of full cycles to avoid aliasing effects
ncoarsegrains = floor(length(AbelData.photonbins)/binspercoarsegrain);
coarsephotons = zeros(1,ncoarsegrains);
for j = 1:ncoarsegrains
    coarsephotons(j) = sum(double(AbelData.photonbins(((j-1)*binspercoarsegrain+1):(j*binspercoarsegrain))));
end
coarsetimes = (0:(ncoarsegrains-1))*binspercoarsegrain*AbelInfo.binlength*1e-9;
figure
plot(coarsetimes,coarsephotons/(coarsetimes(2)));
xlabel('Time (s)')
ylabel('Photon Counts (Hz)')

end
end
end
end