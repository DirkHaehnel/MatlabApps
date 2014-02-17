function AbelData = traceAbelDataFiles(AbelInfo)
%function AbelData = traceAbelDataFiles(AbelInfo)
%Decode the saved data from an ABEL trapping experiment
%Input:
%ABELinfo is the output from importAbelFiles function, which is in general
%the experiimental parameters
%Output:
%photonbins, voltages and startbin as a vector
%AbelData = struct('photonbins',photonbins,'voltages',voltages,'startbin',startbin)
%Dataformat:
%photonbins is a uint8 (Mx1) vector of bin-by-bin photon counts
%voltages is an int16 (Mx2) vector of bin-by-bin voltages - NOT TRUE PHYSICAL VOLTAGES!
%Extra:
%voltages can be converted to physical voltages by dividing by 3276.8
%This step is not performed here in an effort to save file space

fbackbytes = 16;
photonfname = AbelInfo.photonDat;
fbackfname = AbelInfo.feedbackDat;
pathToDatafname=AbelInfo.pathToData;
temp = dir(fullfile(pathToDatafname,fbackfname));
totalfbackbytes = temp.bytes;
clear temp;
photonfile = fopen(fullfile(pathToDatafname,photonfname),'r');
fbackfile = fopen(fullfile(pathToDatafname,fbackfname),'r');    
%Synchronize File Start Positions
[currphotoncycle,currphotonbin,photons] = getnextkalmanphoton(photonfile);
[currcycle,currbin,currest,currvar] = getnextkalmanest_int(fbackfile);
cyclediff = double(currcycle) - double(currphotoncycle);
if cyclediff > 2^15
    cyclediff = cyclediff - 2^16;
elseif cyclediff < -2^15
    cyclediff = cyclediff + 2^16;
end
if cyclediff > 0 || (~cyclediff && currbin >= currphotonbin) %extra photons at the start of the file prior to the first feedback data
    while cyclediff > 0 || (~cyclediff && currbin >= currphotonbin)
        [currphotoncycle,currphotonbin,photons] = getnextkalmanphoton(photonfile);
        cyclediff = double(currcycle) - double(currphotoncycle);
        if cyclediff > 2^15
            cyclediff = cyclediff - 2^16;
        elseif cyclediff < -2^15
            cyclediff = cyclediff + 2^16;
        end
    end

else 
    while cyclediff < 0 || (~cyclediff && currbin < currphotonbin)
        [currcycle,currbin,currest,currvar] = getnextkalmanest_int(fbackfile);

        cyclediff = double(currcycle) - double(currphotoncycle);
        if cyclediff > 2^15
            cyclediff = cyclediff - 2^16;
        elseif cyclediff < -2^15
            cyclediff = cyclediff + 2^16;
        end
    end

    while cyclediff > 0 || (~cyclediff && currbin >= currphotonbin)
        [currphotoncycle,currphotonbin,photons] = getnextkalmanphoton(photonfile);
        cyclediff = double(currcycle) - double(currphotoncycle);
        if cyclediff > 2^15
            cyclediff = cyclediff - 2^16;
        elseif cyclediff < -2^15
            cyclediff = cyclediff + 2^16;
        end
    end

end

%At this point, [currcycle,currbin,currest,currvar] hold the first usable
%feedback point, and [currphotoncycle,currphotonbin,photons] hold the first usable
%photon point

varscale = 4;

dt = AbelInfo.binlength * 1e-9; %to avoid repeated conversion to seconds from ns

inversemobility = int32(round(AbelInfo.EODscale/(AbelInfo.mobility*dt)));
inversediffconst = uint32(round(2^varscale*AbelInfo.beamwaist^2/(2*AbelInfo.diffconst*dt)));
EODs = int32(round([AbelInfo.EODx,AbelInfo.EODy]*3276.8));
maxstep = int32(round(AbelInfo.maxV*3276.8/abs(inversemobility)));

[nextfbackcycle,nextfbackbin,nextest,nextvar] = getnextkalmanest_int(fbackfile);
if currbin - nextfbackbin ~= 1 && ABELinfo.binspercycle + currbin - nextfbackbin ~=1
    fprintf(['Error: feedback file appears to be incorrect.\n'...
        'Current cycle %d, bin %d; next cycle %d, bin %d\n'],[currcycle,currbin,nextfbackcycle,nextfbackbin]);
end

if ~AbelInfo.mobilityon
    fprintf('Error: Mobility must be on!\n')
end

startvoltages = zeros(3,2,'int32');

[currcycle,currbin] = increment(currcycle,currbin,AbelInfo.binspercycle);
startvoltages(1,:) = min(maxstep,max(-maxstep,currest)); %Guess some start startvoltages
startvoltages(2,:) = min(maxstep,max(-maxstep,currest - startvoltages(1,:)));
startvoltages(3,:) = min(maxstep,max(-maxstep,currest - (startvoltages(1,:)+startvoltages(2,:))));
goodcycles = 0;

%We're now synchronized, so we can start recording 
%The first three voltages should already be correct

numbins = (1 + (totalfbackbytes - ftell(fbackfile))/fbackbytes) * (AbelInfo.binspercycle - 1);
%This is an upper bound on the number of bins of the trace, to avoid copying data
startbin = currbin;
fclose(fbackfile);

photonbins = zeros(numbins,1,'uint8');
voltages = zeros(numbins+3,2,'int32');

voltages(1:3,:) = startvoltages;

for i = 1:numbins
    if ~mod(i,100000)
        fprintf('Bin %d of %d\n',i,numbins)
    end
    currest = currest - voltages(i,:);
    currvar = currvar + 2^varscale;
    voltages(i+3,:) = min(maxstep,max(-maxstep,currest - (voltages(i+1,:)+voltages(i+2,:))));
    if currcycle == currphotoncycle && currbin == currphotonbin
        photonbins(i) = photons;
        currest = int32(inversediffconst)*currest + int32(uint32(photons)*currvar)*EODs(currbin,:);
        currest = int32(fix(double(currest)/double(inversediffconst+uint32(photons)*currvar)));
        currvar = uint32(floor(double(currvar*inversediffconst)/double(inversediffconst+uint32(photons)*currvar)));
        [currphotoncycle,currphotonbin,photons] = getnextkalmanphoton(photonfile);
        if feof(photonfile) %stop immediately when reached end of photon file
            break
        end
    end
    [currcycle,currbin] = increment(currcycle,currbin,AbelInfo.binspercycle);
end

fclose(photonfile);

photonbins = photonbins(1:i);
%voltages = voltages(1:i,:)*inversemobility/3276.8;
%uncoment this fr real voltage values and coment the line below
voltages = int16(voltages(1:i,:)*inversemobility);
AbelData = struct('photonbins',photonbins,'voltages',voltages,'startbin',startbin);

