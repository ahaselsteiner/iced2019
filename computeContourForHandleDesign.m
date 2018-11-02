%% computeContourForHandleDesign
% Computes a 99% highest density contour with the variables hand length
% and hand width. This script was used for the work presented at ICED 2019.

%% Dependencies
% * Matlab 2017b with the Statistics and Machine Learning Toolbox.
% * compute-hdc version 1.10 (can be downloaded from 
% https://github.com/ahaselsteiner/compute-hdc/releases/tag/1.1.0).

X1_GRID_RESOLUTION = 0.2; % Reduce the resolution to compute faster.
X2_GRID_RESOLUTION = 0.1; % Reduce the resolution to compute faster.
START_X1 = 100;
END_X1 = 300;
START_X2 = 50;
END_X2 = 150;
LABEL_X1 = 'hand length (mm)';
LABEL_X2 = 'hand breadth (mm)';
BANDWIDTH_FACTOR_SILVERMAN = 1;

female = importAnsurFile('ANSUR_II_FEMALE_Public.csv', 2, 21);
male = importAnsurFile('ANSUR_II_MALE_Public.csv', 2, 21);

combined.subjectId = [female.subjectId; male.subjectId];
combined.handBreadth = [female.handBreadth; male.handBreadth];
combined.handLength = [female.handLength; male.handLength];
combined.stature = [female.stature; male.stature];
combined.weightKg = [female.weightKg; male.weightKg];
combined.gender = [female.gender; male.gender];
combined.date = [female.date; male.date];
combined.subjectsBirthLocation = [female.subjectsBirthLocation; male.subjectsBirthLocation];
combined.age = [female.age; male.age];

% Choose the bandwidth according to Silverman's rule of thumb. See 
% https://de.mathworks.com/help/stats/mvksdensity.html for an explanation.
dataX1 = combined.handLength;
dataX2 = combined.handBreadth;
sigma = std([dataX1 dataX2]);
n = length(dataX1);
d = 2;
b = BANDWIDTH_FACTOR_SILVERMAN * sigma*(4/((d+2)*n))^(1/(d+4));

% --- Density estimation, f_{DE} ---
x1Ticks = START_X1:X1_GRID_RESOLUTION:END_X1;
x2Ticks = START_X2:X2_GRID_RESOLUTION:END_X2;
[x1Grid x2Grid] = meshgrid(x1Ticks, x2Ticks);
x1Vector = x1Grid(:);
x2Vector = x2Grid(:);
xi = [x1Vector x2Vector];
[f xi bw] = mvksdensity([dataX1 dataX2], xi, 'bandwidth', b);
[F xi bw] = mvksdensity([dataX1 dataX2], xi, 'bandwidth', b, 'function', 'cdf');

counter = 1;
for i = 1:length(x1Ticks)
    for j = 1:length(x2Ticks)
        fGrid(j,i) = f(counter);
        FGrid(j,i) = F(counter);
        counter = counter + 1;
    end
end

dX1 = X1_GRID_RESOLUTION;
dX2 = X2_GRID_RESOLUTION;
gridCenterPoints = {(x1Grid(1,1) + 0.5*dX1) : dX1 : x1Grid(1,end); ...
    (x2Grid(1,1) + 0.5*dX2) : dX2 : x2Grid(end,1)};

% Build a model using the required data structure.
M.name = 'ANSUR dataset with KDE';
M.modelType = 'KDE';
M.cdf= FGrid';
M.pdf = fGrid';
M.cdfGrid = {x1Grid'; x2Grid'};
M.gridCenterPoints = gridCenterPoints;
M.labels = {LABEL_X1; LABEL_X2};

% --- Contour construction, f_{CC} ---
[fm, x1Hdc, x2Hdc] = computeHdc(M, 0.01, M.gridCenterPoints, 0);

% --- Design condition selection, f_{DC} ---
[value i] = max(x1Hdc{1});
indices = i;
[value i] = max(x2Hdc{1});
indices = [indices i];
[value i] = min(x1Hdc{1});
indices = [indices i];
[value i] = min(x2Hdc{1});
indices = [indices i];
[value i] = max((x1Hdc{1} / mean(x1Hdc{1})) .* (x2Hdc{1} / mean(x2Hdc{1})));
indices = [indices i];
[value i] = min((x1Hdc{1}  / mean(x1Hdc{1})) .* (x2Hdc{1}  / mean(x2Hdc{1})));
indices = [indices i];

% Select the design conditions along the design contour.
x1DC = x1Hdc{1}(indices);
x2DC = x2Hdc{1}(indices);

% Plot the sample and the design contour
fig = figure();
plot(female.handLength, female.handBreadth, '+k')
hold on
plot(male.handLength, male.handBreadth, 'ok')
plot(x1Hdc{1}, x2Hdc{1}, '-b')
plot(x1DC, x2DC, 'xk', 'MarkerSize', 12, 'LineWidth', 2)
xlabel(LABEL_X1)
ylabel(LABEL_X2)
xlim([140 250])
ylim([55 110])
legend('female', 'male', '99% highest density contour', 'design conditon', 'location', 'southeast');
legend boxoff;
set(fig, 'Units', 'centimeters', 'Position', [0, 0, 10, 8], 'PaperUnits', 'Inches', 'PaperSize', [10 8])

