function AbelDataTrapped = findTrappedMolecules(photonbins,binspercycle,binlength)
%function AbelDataTrapped = findTrappedMolecules(photonbins,binspercycle,binlength)
%The input comes from the traceAbelData function
%Input:
%photonbins, binspercycle and binlength of that oparticular
%datafile
%Important: Must be from same file otherwise stupid idea!
%Output:
%AbelDataTrapped = struct('trapstarts',trapstarts,'traplengths',traplengths);
%A struct of the position where a trapping starts and the length of that
%trapping event
%AbelDataTrapped returns a struct with trapstart and traplength is the resultvector of the maximum likelood analysis

niter = 100; 
%probably excessive, but it goes so fast it doesn't really matter

dt = binlength*1e-9;
%Sum over full cycles to avoid fluctuations due to laser motion
cycletime = binspercycle*dt;
ncycles = floor(length(photonbins)/binspercycle);
photoncycles = photonbins(1:binspercycle:binspercycle*ncycles);
for i = 2:binspercycle;
    photoncycles = photoncycles + photonbins(i:binspercycle:binspercycle*ncycles);
end

countrate0 = zeros(1,niter+1);
countrate1 = zeros(1,niter+1);
k01 = zeros(1,niter+1);
k10 = zeros(1,niter+1);
countrate0(1) = 6e3*cycletime; % Starting value for the background count rate
countrate1(1) = 4e4*cycletime; % Starting value for the signal+background count rate
k01(1) = 0.1*cycletime; % Starting value for the probability for a molecule to enter
k10(1) = 0.1*cycletime; % Starting value for the probability for a molecule to exit
for i = 1:niter
    %call  the compiled C_twostateBW c code function
    %run mex C_twostateBW.c in C:\tmp since compiler doesn't liek space's
    [countrate0(i+1),countrate1(i+1),k01(i+1),k10(i+1)] = C_twostateBW(countrate0(i),countrate1(i),k01(i),k10(i),double(photoncycles));
end
 %call  the compiled C_twostateV c code function
    %run mex C_twostateV.c in C:\tmp since compiler doesn't liek space's
[LL,states] = C_twostateV(countrate0(niter+1),countrate1(niter+1),k01(niter+1),k10(niter+1),double(photoncycles));

trapstarts = 1 + binspercycle*(find(diff([0;states])==1) - 1);
% Indices in photonbins when a molecule has first entered
traplengths = 1 + binspercycle*(find(diff([states;0])==-1)) - trapstarts;
% Length of time in bins (not cycles) that the molecule remains trapped
% So the first molecule's counts would be photonbins(trapstarts(1):trapstarts(1)+traplengths(1)-1)
AbelDataTrapped = struct('trapstarts',trapstarts,'traplengths',traplengths);
