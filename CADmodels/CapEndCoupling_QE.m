%% Maps for the external quality factor
clear
clc

PAF = PathsAndFiles;

AES.f0 = 1.575;
AES.threshold = 0.001;
AES.MaxFunEvals = 20;

hp = 6:4:22;
sn = 1;

[X1,X2] = ndgrid(hp,sn);

hpi = reshape(X1, [], 1);
sni = reshape(X2, [], 1);

nCases = numel(X1);

% Optimization settings
x1 = 0.40; x2 = 0.60;
CheckTuning = @(x, optimValues, state) outfunGeneral(x, optimValues, state, AES);
f2opt = @(x) tuning(x,PAF,AES);
options = optimset("MaxFunEvals",AES.MaxFunEvals,'OutputFcn',CheckTuning);

Wtn1 = zeros(size(hpi));
QE = zeros(size(hpi));
fresload = zeros(size(hpi));

for idx = 1:nCases

    mws = PAF.mws;
    mws.invoke('StoreDoubleParameter','hp', hpi(idx));
    mws.invoke('StoreDoubleParameter','sn',  sni(idx));
    mws.invoke('Rebuild'); % Rebuild a structure
    mws.invoke('save');
    clear mws

    [Wtn1(idx),fval,exitflag,output] = fminbnd(f2opt,x1,x2,options);

    Out = OutputSim(Wtn1(idx),PAF);
    fresload(idx) = Out.fresload;
    QE(idx) = Out.Qext;

    save("CapEndCoupling_QE.mat", "Wtn1", "QE", "fresload")

end

invoke(PAF.cst, 'quit');

%% Functions

function stop = outfunGeneral(x, optimValues, state, AES)
stop = false;

if optimValues.fval <= AES.threshold
    stop = true;
end
end

function PAF = PathsAndFiles

PAF.folder = "C:\Users\giannetti\OneDrive - unifi.it\Università\Dottorato\Lavori\FiltroE1Galileo\MATLAB\CST_API\";
PAF.filename = "CapEndCoupling_QE.cst";
PAF.file2open = strcat(PAF.folder, PAF.filename);

PAF.cst = actxserver('CSTStudio.Application');
PAF.mws = invoke(PAF.cst , 'NewMWS');
PAF.mws.invoke('OpenFile',PAF.file2open);

end

function cost = tuning(x,PAF,AES)

fprintf("Normalized screw penetration: %.4f\n", x)

Out = OutputSim(x,PAF);

cost = abs(Out.fresload - AES.f0);

fprintf("\n")

end

function Out = OutputSim(x,PAF)

mws = PAF.mws;

mws.invoke('StoreDoubleParameter','Wtn1', x); % Send into parameter theta a new value
mws.invoke('Rebuild'); % Rebuild a structure

mws.invoke('save');

hSolver = mws.invoke('EigenmodeSolver');
checkSim = invoke(hSolver,'Start');
if ~checkSim
    warning("The simulation failed")
end

Nmodes   = invoke(hSolver,'GetNumberOfModesCalculated');
fresload = invoke(hSolver,'GetLoadedFrequencyInHz',1)*1e-9;
Qext     = invoke(hSolver,'GetModeExternalQFactor',1);

mws.invoke('save');

clear mws hSolver

Out.Nmodes = Nmodes;
Out.fresload = fresload;
Out.Qext = Qext;

fprintf("Loaded resonant frequency: %.3f\n", fresload)
fprintf("External quality factor: %.2f\n", Qext)

end