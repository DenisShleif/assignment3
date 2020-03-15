%% Main Code

function [] = assignment3Code_DenisShleifman_101001778()
clear
clc
close all

const.me = 9.10938356e-31; %electron mass(kg)
const.re = 2.8179403227e-15; %electron radius (m)
const.kb = 1.38066e-23; %Boltzmans Constant (J/K)
const.T = 300; %Board Temperature (K)
const.tmn = 0.2e-12; % Mean time between collisions (s)
const.e = -1.60217662e-19; % elemenrary charge (C)
const.n = 1e19; %Particle denistiy (m^-2) 

sim.particles = 1000; %number of particles to initialize in the simulation
sim.delay = 0.01; %time in seconds between executing steps
sim.steps = 1000; %number of iterations to perform
sim.movie = false; %set to true if you would like it to loop through each iteration
sim.scatter = true; %set whether boundaries cause
sim.heatmapBins = [10,10]; %number of bins for 2D heatmap [x,y]
sim.bottleNeck = false; %set to true to place the bottleneck
sim.solveEField = false; %set to true to solve E field

board.xMin = 0; %Minimum X coordinate for the board (m)
board.yMin = 0; %Minimum Y coordinate for the board (m)
board.xMax = 200e-9; %Maximum X coordinate for the board (m)
board.yMax = 100e-9; %Maximum Y coordinate for the board (m)
board.Lb = 50e-9; %Length of the Bottle Neck (nm)
board.Wb = 10e-9; %With of the Bottle Neck (nm)
board.sigBox = 1e-2; %Conductivity of forbidden region (S)
board.sigOut = 1; %Conductivity of non forbidden region (S)
board.h = 5e-9; %Grid for solving for Electric Field
board.borderV = [0.1,0,0,0]; %Voltage of Borders [Left,Right,Bottom,Top]
board.borderFloat = [false,false,true,true]; %Set to true to float certain borders [Left,Right,Bottom,Top]

part.colour = blanks(sim.particles);
part.colour(1) = 'b';
part.colour(2) = 'g';
part.colour(3) = 'r';
part.colour(4) = 'c';
part.colour(5) = 'm';
part.colour(6) = 'y';
part.colour(7) = 'k';

sim.dist.mult = 1e9;
sim.dist.unit = 'nm';
sim.time.mult = 1e12;
sim.time.unit = 'ps';
sim.speed.mult = 1e-3;
sim.speed.unit = 'km/s';
sim.curr.mult = 1e3;
sim.curr.unit = 'mA';

performSimuation(const,part,board,sim);
end

function [] = performSimuation(const,part,board,sim)
[const,part,board,sim,stat,plt] = initialize(const,part,board,sim);

for n = 1:sim.steps
    [part,stat,plt,sim] = performIteration(const,part,board,sim,stat,plt,n);
end
updatePlot(plt,part,stat,sim,sim.steps,true);
outputCalcResults(const,sim,stat,board);
plotHeatMaps(board,sim,stat);
end

%% Initialization of all variables

function [const,part,board,sim,stat,plt] = initialize(const,part,board,sim)
const = initializeConst(const);
board = initializeBoard(board,const,sim);
sim = initializeSim(sim,const,board.minSideLength);
part = initializePart(part,board,sim,const);
stat = initializeStats(board,sim);
plt = setupPlot(board,sim,const);
end

function [const] = initializeConst(const)
const.meff = const.me*0.26;
const.de = const.re*2;
const.vth = sqrt((2*const.T*const.kb)/const.meff);
const.vav = const.vth*sqrt(pi/4);
const.MFP = const.vav *const.tmn;
end

function [sim] = initializeSim(sim,const,minSideLength)
sim.dt = minSideLength/(100*const.vth);
sim.minStep = minSideLength/100;
sim.outputPercent = 0;
sim.Pscat = 1 - exp(-sim.dt/const.tmn);
sim.cornerPrec = const.vth*sim.dt;
end

function [part] = initializePart(part,board,sim,const)
part = generateInitialPosition(part,board,sim);
part = assignVelocities(part,sim,const);
part.colTime = zeros(1,sim.particles);
part.colDist = zeros(1,sim.particles);
part.v = zeros(1,sim.particles);
part.v2 = zeros(1,sim.particles);
part.temp = zeros(1,sim.particles);
part.traceIndex = part.colour ~= ' ';
end

function [part] = assignVelocities(part,sim,const,I)
if nargin == 3
    I = true(1,sim.particles);
end

part.vx(I) = nrmrnd(0,1,1,sum(I))*const.vth/sqrt(2);
part.vy(I) = nrmrnd(0,1,1,sum(I))*const.vth/sqrt(2);
end

function [out] = generatePoints(N,Min,Max)
out = rand(1,N)*(Max - Min) + Min;
end

function nrmmatrix = nrmrnd(mu, sigma, sz1, sz2)
nrmmatrix = mu+sigma*randn(sz1,sz2);
end

function [part] = generateInitialPosition(part,board,sim)
I = true(1,sim.particles);
while sum(I) ~= 0
    part.x(I) = generatePoints(sum(I),board.xMin,board.xMax);
    part.y(I) = generatePoints(sum(I),board.yMin,board.yMax);
    I = false(1,sim.particles);
    if sim.bottleNeck
        for n = 1:length(board.boxes)
            I = I | isInsideBox(board.boxes{n},part.x,part.y);
        end
    end
end
end

function [I] = isInsideBox(box,x,y)
I = (x > box(1) & x < box(3) & y > box(2) & y < box(4));
end

function [stat] = initializeStats(board,sim)
stat.time = zeros(1,sim.steps);
stat.temp = zeros(1,sim.steps);
stat.avtemp = zeros(1,sim.steps);
stat.MFP = zeros(1,sim.steps);
stat.tmn = zeros(1,sim.steps);
stat.driftCurr = zeros(1,sim.steps);
stat.avDriftCurr = zeros(1,sim.steps);
stat.avV = zeros(1,sim.steps);
stat.avVrms = zeros(1,sim.steps);
stat.heatMap = generateHeatMapXY(board,sim);
stat.avHeatMap = zeros(sim.heatmapBins(1)+1,sim.heatmapBins(2)+1);
stat.avTempMap = zeros(sim.heatmapBins(1)+1,sim.heatmapBins(2)+1);
end

function [heatMap] = generateHeatMapXY(board,sim)
heatMap.x = linspace(board.xMin,board.xMax,sim.heatmapBins(1)+1);
heatMap.y = linspace(board.yMin,board.yMax,sim.heatmapBins(2)+1);
[heatMap.X,heatMap.Y] = meshgrid(heatMap.x*sim.dist.mult,heatMap.y*sim.dist.mult);
end
%% Setup Board


function [board] = initializeBoard(board,const,sim)
board.Lx = board.xMax - board.xMin;
board.Ly = board.yMax - board.yMin;
board.A = board.Lx*board.Ly;
board.minSideLength = min([board.Lx,board.Ly]);
board = initializeBoundaries(board,sim);
if sim.solveEField
    [board] = initializeEFieldBoard(board);
    [board] = solveEField(board,sim,const);
end
board.dVx = board.borderV(2) - board.borderV(1);
board.dVy = board.borderV(4) - board.borderV(3);
board.Ex = -board.dVx/board.Lx;
board.Ey = -board.dVy/board.Ly;
board.Fx = const.e*board.Ex;
board.Fy = const.e*board.Ey;
board.ax = board.Fx/const.meff;
board.ay = board.Fy/const.meff;
end

function [board] = initializeEFieldBoard(board)
board.x = board.xMin:board.h:board.xMax;
board.Nx = length(board.x);
board.y = board.yMin:board.h:board.yMax;
board.Ny = length(board.y);
[board.X,board.Y] = meshgrid(board.x,board.y);
board.N = board.Nx*board.Ny;
board.sig = zeros(board.Ny,board.Nx);
board.xlim = [(board.Lx - board.Lb)/2,(board.Lx + board.Lb)/2];
board.ylim = [(board.Ly - board.Wb)/2,(board.Ly + board.Wb)/2];
I = board.X >= board.xlim(1) & board.X <= board.xlim(2) & ~(board.Y > board.ylim(1) & board.Y < board.ylim(2));
board.sig(I) = board.sigBox;
board.sig(~I) = board.sigOut;
board.sigAv = mean(board.sig,'all');
end

function [eq] = getEqSpace(nx,yPos,xPos)
eq = xPos + (yPos-1)*nx;
end

function [yPos,xPos] = getPosSpace(nx,eq)
xPos = mod(eq-1,nx)+1;
yPos = floor((eq-1)/nx) + 1;
end

function [G,val] = assignG(board,curr,G,delta,val)
eq = getEqSpace(board.Nx,curr.y + delta(1),curr.x + delta(2));
if eq <= 0 || eq > board.N
    return;
end
if nargin == 4
    val = mean(board.sig(curr.y + [0,delta(1)],curr.x + [0,delta(2)]),'all');
end
G(curr.n,eq) = val;
end

function [G,B] = assignEdge(board,G,B,curr,borderNum)
if board.borderFloat(borderNum)
    counter = 0;
    val = zeros(1,4);
    if curr.x - 1 >= 1
        counter = counter + 1;
        [G,val(counter)] = assignG(board,curr,G,[0,-1]);
    end
    if curr.x + 1 <= board.Nx
        counter = counter + 1;
        [G,val(counter)] = assignG(board,curr,G,[0,+1]);
    end
    if curr.y - 1 >= 1
        counter = counter + 1;
        [G,val(counter)] = assignG(board,curr,G,[-1,0]);
    end
    if curr.y + 1 <= board.Ny
        counter = counter + 1;
        [G,val(counter)] = assignG(board,curr,G,[+1,0]);
    end
    G = assignG(board,curr,G,[0,0],-sum(val(1:counter)));
    B(curr.n) = 0;
else
    G = assignG(board,curr,G,[0,0],1);
    B(curr.n) = board.borderV(borderNum);
end
end

function [board] = solveEField(board,sim,const)
G = spalloc(board.N,board.N,board.N*5);
B = zeros(board.N,1);
val = zeros(1,4);
for n = 1:board.N
    curr.n = n;
    [curr.y,curr.x] = getPosSpace(board.Nx,curr.n);
    if curr.x == 1
        [G,B] = assignEdge(board,G,B,curr,1);
    elseif curr.x == board.Nx
        [G,B] = assignEdge(board,G,B,curr,2);
    elseif curr.y == 1
        [G,B] = assignEdge(board,G,B,curr,3);
    elseif curr.y == board.Ny
        [G,B] = assignEdge(board,G,B,curr,4);
    else
        [G,val(1)] = assignG(board,curr,G,[-1,0]);
        [G,val(2)] = assignG(board,curr,G,[1,0]);
        [G,val(3)] = assignG(board,curr,G,[0,-1]);
        [G,val(4)] = assignG(board,curr,G,[0,+1]);
        G = assignG(board,curr,G,[0,0],-sum(val));
        B(curr.n) = 0;
    end
end

V = mldivide(G,B);
V = reshape(V,[board.Nx,board.Ny])';
[Ex,Ey] = gradient(V,board.h);
Ex = -Ex;
Ey = -Ey;
Ex(isnan(Ex)) = 0;
Ey(isnan(Ey)) = 0;
Fx = const.e*Ex;
Fy = const.e*Ey;
Ax = Fx/const.meff;
Ay = Fy/const.meff;

board.Ax = griddedInterpolant(board.X',board.Y',Ax','spline');
board.Ay = griddedInterpolant(board.X',board.Y',Ay','spline');

figure;
surf(board.X.*sim.dist.mult,board.Y.*sim.dist.mult,V);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Voltage (V)');
title('Voltage Plot');
grid on;
view(25,28)

board.xMin = 0; %Minimum X coordinate for the board (m)
board.yMin = 0; %Minimum Y coordinate for the board (m)
board.xMax = 200e-9; %Maximum X coordinate for the board (m)
board.yMax = 100e-9; %Maximum Y coordinate for the board (m)

figure;
ax = subplot(1,1,1);
quiver(board.X.*sim.dist.mult,board.Y.*sim.dist.mult,Ex,Ey,0.4);
hold on
drawBoundaries(ax,sim,board)
axis([board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
xlabel('x (nm)');
ylabel('y (nm)');
title('Electric Field Plot');
grid on

figure;
ax = subplot(1,1,1);
quiver(board.X'.*sim.dist.mult,board.Y'.*sim.dist.mult,board.Ax(board.X',board.Y'),board.Ay(board.X',board.Y'),0.4);
hold on
drawBoundaries(ax,sim,board)
axis([board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
xlabel('x (nm)');
ylabel('y (nm)');
title('Particle Acceleration Plot');
grid on
end

function [board] = initializeBoundaries(board,sim)
board.startX = board.xMin - board.Lx/5;
board.endX = board.xMax + board.Lx/5;

if ~sim.bottleNeck
    board.boundaries{1} = [board.startX,board.yMin,board.endX,board.yMin,1];
    board.boundaries{2} = [board.startX,board.yMax,board.endX,board.yMax,-1];
end

board.xBoxLeft = board.xMin + (board.Lx - board.Lb)/2;
board.xBoxRight = board.xBoxLeft + board.Lb;
board.yBoxTop = board.yMin + (board.Ly - board.Wb)/2;
board.yBoxBot = board.yBoxTop + board.Wb;

board.boxes{1} = [board.xBoxLeft,board.yMin,board.xBoxRight,board.yBoxTop];
board.boxes{2} = [board.xBoxLeft,board.yBoxBot,board.xBoxRight,board.yMax];
board.boundaries{1} = [board.startX,board.yMin,board.xBoxLeft,board.yMin,1];
board.boundaries{2} = [board.xBoxLeft,board.yMin,board.xBoxLeft,board.yBoxTop,-1];
board.boundaries{3} = [board.xBoxLeft,board.yBoxTop,board.xBoxRight,board.yBoxTop,1];
board.boundaries{4} = [board.xBoxRight,board.yBoxTop,board.xBoxRight,board.yMin,1];
board.boundaries{5} = [board.xBoxRight,board.yMin,board.endX,board.yMin,1];
board.boundaries{6} = [board.startX,board.yMax,board.xBoxLeft,board.yMax,-1];
board.boundaries{7} = [board.xBoxLeft,board.yMax,board.xBoxLeft,board.yBoxBot,-1];
board.boundaries{8} = [board.xBoxLeft,board.yBoxBot,board.xBoxRight,board.yBoxBot,-1];
board.boundaries{9} = [board.xBoxRight,board.yBoxBot,board.xBoxRight,board.yMax,1];
board.boundaries{10} = [board.xBoxRight,board.yMax,board.endX,board.yMax,-1];
board.corners{1} = [board.xBoxLeft,board.yMin,2];
board.corners{2} = [board.xBoxRight,board.yMin,1];
board.corners{3} = [board.xBoxLeft,board.yMax,3];
board.corners{4} = [board.xBoxRight,board.yMax,4];
board.corners{5} = [board.xBoxLeft,board.yBoxTop,2];
board.corners{6} = [board.xBoxRight,board.yBoxTop,1];
board.corners{7} = [board.xBoxLeft,board.yBoxBot,3];
board.corners{8} = [board.xBoxRight,board.yBoxBot,4];

for n = 1:length(board.boxes)
    board.A = board.A - abs((board.boxes{n}(3) - board.boxes{n}(1))*(board.boxes{n}(4) - board.boxes{n}(2)));
end
end

%% Setup Plots

function [plt] = setupPlot(board,sim,const)
plt.fig1 = setupParticlePlot(board,sim);
plt.fig2 = setupTempPlot(sim);
plt.fig3 = setupVelocityDist(sim,const);
plt.fig4 = setupCurrDenistyPlot(sim);
end

function [fig] = setupParticlePlot(board,sim)
fig.h = figure;
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
axis(fig.ax1,[board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
hold(fig.ax1, 'on');
if sim.movie
    fig.plotPos = plot(fig.ax1,1,1,'.b');
    fig.plotVel = quiver(fig.ax1,[],[],[],[],2,'b');
end
drawBoundaries(fig.ax1,sim,board);
xlabel(fig.ax1,sprintf('x (%s)',sim.dist.unit));
ylabel(fig.ax1,sprintf('y (%s)',sim.dist.unit));
title(fig.ax1,sprintf('Particle Trajectories (7/%d)',sim.particles));
grid(fig.ax1,'on');
end

function [] = drawBoundaries(ax,sim,board)
if sim.bottleNeck
    for n = 1:length(board.boundaries)
        currentBoundary = board.boundaries{n};
        line(ax,[currentBoundary(1),currentBoundary(3)]*sim.dist.mult,[currentBoundary(2),currentBoundary(4)]*sim.dist.mult,'Color','k');
    end
end
end

function [fig] = setupTempPlot(sim)
fig.h = figure();
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
hold(fig.ax1, 'on')
fig.plotTemp = plot(fig.ax1,1,1,'-b');
fig.plotAvTemp = plot(fig.ax1,1,1,'-r');
legend(fig.ax1,{'Current','Average'});
xlabel(fig.ax1,sprintf('Time (%s)',sim.time.unit));
ylabel(fig.ax1,'Temperature (K)');
grid(fig.ax1,'on');
end

function [fig] = setupCurrDenistyPlot(sim)
fig.h = figure();
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
hold(fig.ax1, 'on')
fig.plotCurr = plot(fig.ax1,1,1,'-b');
fig.plotAv = plot(fig.ax1,1,1,'-r');
legend(fig.ax1,{'Current','Average'});
xlabel(fig.ax1,sprintf('Time (%s)',sim.time.unit));
ylabel(fig.ax1,sprintf('Current (%s)',sim.curr.unit));
grid(fig.ax1,'on');
end

function [fig] = setupVelocityDist(sim,const)
fig.h = figure();
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
fig.axis = axis(fig.ax1);
hold(fig.ax1, 'on')
fig.vHist = histogram(fig.ax1,[],'FaceColor','g');
xline(fig.ax1,const.vth*sim.speed.mult,'Label','vth','Color','r','LineWidth',2);
fig.Vav = xline(fig.ax1,0,'Label','Average Velocity','Color','b','LineWidth',2,'LabelVerticalAlignment','bottom');
fig.Vrms = xline(fig.ax1,0,'Label','Vrms','Color','k','LineWidth',2,'LabelVerticalAlignment','middle');
xlabel(fig.ax1,sprintf('Velocity (%s)',sim.speed.unit));
ylabel(fig.ax1,'Counts');
grid(fig.ax1,'on');
end

%%Perform Iteration

function [part,stat,plt,sim] = performIteration(const,part,board,sim,stat,plt,n)
lastX = part.x;
lastY = part.y;

part = boundaryCondition(part,board,sim,const);
part = calcValues(part,const,sim.dt);
part = updateAdd(part,lastX,lastY,board.Lx);
part = scatter(part,sim,const);
stat = updateStats(stat,const,board,part,sim.dt,n);
plt = updatePlot(plt,part,stat,sim,n);
if sim.movie
    pause(sim.delay);
else
    sim = outputPercentDone(sim,n);
end
end

%% Boundary Conditions

function [part] = boundaryCondition(part,board,sim,const)
dt = ones(1,sim.particles) * sim.dt;
forbid = true(1,sim.particles);
counter = 0;
if sim.solveEField
    part.ax = board.Ax(part.x,part.y);
    part.ay = board.Ay(part.x,part.y);
else
    part.ax = ones(1,sim.particles)*board.ax;
    part.ay = ones(1,sim.particles)*board.ay;
end
while sum(forbid) > 0
    [x,y,vx,vy] = calcNewPos(part,dt);
    forbid = checkForbiddenRedions(sim,board,x,y);
    part.x(~forbid) = x(~forbid);
    part.y(~forbid) = y(~forbid);
    part.vx(~forbid) = vx(~forbid);
    part.vy(~forbid) = vy(~forbid);
    dt(~forbid) = 0;
    if sim.bottleNeck
        for n = 1:length(board.corners)
            part = checkCorners(part,sim,const,forbid,board.corners{n});
        end
    end
    for n = 1:length(board.boundaries)
        [part,dt] = lineIntersection(board.boundaries{n},part,sim,const,dt,forbid,x,y);
    end
    counter = counter + 1;
    if counter == 50
        part = reassign(board,part,sim,forbid);
        break;
    end
end
part.x = wrap(part.x,board.xMin,board.xMax);
end

function [sim] = outputPercentDone(sim,currStep)
while (currStep/sim.steps)*100 > sim.outputPercent
    sim.outputPercent = sim.outputPercent + 1;
    fprintf('Percent Complete: %d %%\n',sim.outputPercent);
end
end

function [part] = reassign(board,part,sim,forbid)
buffer = sim.minStep/1000;
moveStep = 2*sim.minStep;

I = (part.x <= (board.xBoxLeft + buffer)) & (part.x >= (board.xBoxLeft - buffer)) & ((part.y <= (board.yBoxTop + buffer)) | (part.y >= (board.yBoxBot - buffer))) & forbid;
part.x(I) = part.x(I) - moveStep;
I = (part.x <= (board.xBoxRight + buffer)) & (part.x >= (board.xBoxRight - buffer)) & ((part.y <= (board.yBoxTop + buffer)) | (part.y >= (board.yBoxBot - buffer))) & forbid;
part.x(I) = part.x(I) + moveStep;

I = (part.x >= (board.xBoxLeft - buffer)) & (part.x <= (board.xBoxRight + buffer)) & (part.y >= (board.yBoxTop - buffer)) & (part.y <= (board.yBoxTop + buffer)) & forbid;
part.y(I) = part.y(I) + moveStep;
I = (part.x >= (board.xBoxLeft - buffer)) & (part.x <= (board.xBoxRight + buffer)) & (part.y >= (board.yBoxBot - buffer)) & (part.y <= (board.yBoxBot + buffer)) & forbid;
part.y(I) = part.y(I) - moveStep;

I = (part.y <= (board.yMin + buffer)) & forbid;
part.y(I) = part.y(I) + moveStep;
I = (part.y >= (board.yMin - buffer)) & forbid;
part.y(I) = part.y(I) - moveStep;   
%fprintf('(%f,%f)\n',[part.x(forbid);part.y(forbid)]*sim.dist.mult);
end

function [x,y,vx,vy] = calcNewPos(part,dt)
dvx = part.ax.*dt;
dvy = part.ay.*dt;
vx = part.vx + dvx;
vy = part.vy + dvy;
dx = part.vx.*dt + 0.5*part.ax.*dt.^2;
dy = part.vy.*dt + 0.5*part.ay.*dt.^2;
x = part.x + dx;
y = part.y + dy;
end

function [forbid,x,y] = checkForbiddenRedions(sim,board,x,y,checkBottleNeck)
if nargin == 4
    checkBottleNeck = sim.bottleNeck;
end

forbid = y > board.yMax | y < board.yMin;
if checkBottleNeck
    for n = 1:length(board.boxes)
        forbid = forbid | isInsideBox(board.boxes{n},x,y);
    end
end
end

function [part] = checkCorners(part,sim,const,forbid,currentCorner)
intersect(forbid) = (part.x(forbid) - currentCorner(1)).^2 + (part.y(forbid) - currentCorner(2)).^2 < sim.cornerPrec.^2;
if sum(intersect) == 0
    return;
end
part.x(intersect) = currentCorner(1);
part.y(intersect) = currentCorner(2);
[xSign,ySign] = getExpectedSign(currentCorner);

if sim.scatter
    part = assignScatterBoundaryVelCorner(part,sim,const,intersect,xSign,ySign);
else
    part.vx(intersect) = xSign*abs(part.vy(intersect));
    part.vy(intersect) = ySign*abs(part.vx(intersect));
end
end

function [xSign,ySign] = getExpectedSign(currentCorner)
switch currentCorner(3)
    case 1
        xSign = 1;
        ySign = 1;
    case 2
        xSign = -1;
        ySign = 1;
    case 3
        xSign = -1;
        ySign = -1;
    case 4
        xSign = 1;
        ySign = -1;
end
end

function [part] = assignScatterBoundaryVelCorner(part,sim,const,intersect,xSign,ySign)
IAssign = intersect;
while sum(IAssign) ~= 0
    part = assignVelocities(part,sim,const,IAssign);
    IAssign(intersect) = sign(part.vx(intersect)) ~= xSign & sign(part.vy(intersect)) ~= ySign;
end
end

function [part,dt] = lineIntersection(currentBoundary,part,sim,const,dt,forbid,x,y)
if currentBoundary(1) == currentBoundary(3) %Vertical
    vert = true;
    I = min([part.x;x]) <= currentBoundary(1) & max([part.x;x]) >= currentBoundary(1) & forbid;
    tmp1 = -part.vx(I)./part.ax(I);
    tmp2 = sqrt(part.vx(I).^2 + 2.*part.ax(I).*(currentBoundary(1) - part.x(I)))./part.ax(I);
elseif currentBoundary(2) == currentBoundary(4) %Horizontal
    vert = false;
    I = min([part.y;y]) <= currentBoundary(2) & max([part.y;y]) >= currentBoundary(2) & forbid;
    tmp1 = -part.vy(I)./part.ay(I);
    tmp2 = sqrt(part.vy(I).^2 + 2.*part.ay(I).*(currentBoundary(2) - part.y(I)))./part.ay(I);
end
if sum(I) == 0
    return;
end
t1 = tmp1 + tmp2;
t2 = tmp1 - tmp2;
dtLoss = zeros(1,length(tmp1));
It1 = t1 >= 0 & t1 <= sim.dt & imag(t1) == 0;
dtLoss(It1) = t1(It1);
It2 = t2 >= 0 & t2 <= sim.dt & imag(t2) == 0;
dtLoss(It2) = t2(It2);
dt(I) = dt(I) - dtLoss;
dvx = part.ax(I).*dtLoss;
dvy = part.ay(I).*dtLoss;
part.vx(I) = part.vx(I) + dvx;
part.vy(I) = part.vy(I) + dvy;
if vert
    part.x(I) = currentBoundary(1);
    dy = part.vy(I).*dtLoss + 0.5*part.ay(I).*dtLoss.^2;
    part.y(I) = part.y(I) + dy;
    if sim.scatter
        part = assignScatterBoundaryVel(part,sim,const,I,1,currentBoundary(5));
    else
        part.vx(I) = -part.vx(I);
    end
else
    part.y(I) = currentBoundary(2);
    dx = part.vx(I).*dtLoss + 0.5*part.ax(I).*dtLoss.^2;
    part.x(I) = part.x(I) + dx;
    if sim.scatter
        part = assignScatterBoundaryVel(part,sim,const,I,0,currentBoundary(5));
    else
        part.vy(I) = -part.vy(I);
    end
end
end

function [part] = assignScatterBoundaryVel(part,sim,const,intersect,orientation,expectedSign)
IAssign = intersect;
while sum(IAssign) ~= 0
    part = assignVelocities(part,sim,const,IAssign);
    if orientation == 1
        IAssign(intersect) = sign(part.vx(intersect)) ~= expectedSign;
    elseif orientation == 0
        IAssign(intersect) = sign(part.vy(intersect)) ~= expectedSign;
    end
end
end

function [val] = wrap(val,Min,Max)
val = mod(val,Max - Min) + Min;
end

%% Perform Updates

function [part] = calcValues(part,const,dt)
part.colTime = part.colTime + dt;
part.v2 = part.vx.^2 + part.vy.^2;
part.v = sqrt(part.v2);
part.colDist = part.colDist + dt * part.v;
part.temp = part.v2.*const.meff/(2*const.kb);
end

function [part] = updateAdd(part,lastX,lastY,lengthX)
indexToPlot = part.traceIndex & (abs(lastX - part.x) < (lengthX/2) );
part.addX = [lastX(indexToPlot);part.x(indexToPlot)];
part.addY = [lastY(indexToPlot);part.y(indexToPlot)];
part.addC = part.colour(indexToPlot);
end

function [part] = scatter(part,sim,const)
I = (rand(1,sim.particles) < sim.Pscat);
part.colDist(I) = 0;
part.colTime(I) = 0;
part = assignVelocities(part,sim,const,I);
end

function [stat] = updateStats(stat,const,board,part,dt,currentStep)
stat.time(currentStep) = dt*(currentStep-1);
stat.temp(currentStep) = mean(part.temp);
stat.avtemp(currentStep) = mean(stat.temp(1:currentStep));
stat.MFP(currentStep) = mean(part.colDist);
stat.tmn(currentStep) = mean(part.colTime);
stat.driftCurr(currentStep) = mean(part.vx*const.e*const.n*board.Lx);
stat.avDriftCurr(currentStep) = mean(stat.driftCurr(1:currentStep));
stat.avV(currentStep) = mean(part.v);
stat.avVrms(currentStep) = sqrt(mean(part.v2));
[count,temp] = generateTemperatureAverage(stat,part);
stat.avHeatMap = stat.avHeatMap + count;
stat.avTempMap = stat.avTempMap + temp;
end

function [count,temp] = generateTemperatureAverage(stat,part)
count = zeros(length(stat.heatMap.x),length(stat.heatMap.y));
temp = zeros(length(stat.heatMap.x),length(stat.heatMap.y));
for n = 1:(length(stat.heatMap.x)-1)
    for m = 1:(length(stat.heatMap.y)-1)
        I = part.x > stat.heatMap.x(n) & part.x < stat.heatMap.x(n+1) & part.y > stat.heatMap.y(m) & part.y < stat.heatMap.y(m+1);
        count(n,m) = sum(I);
        if count(n,m) == 0
            temp(n,m) = 0;
        else
            temp(n,m) = mean(part.temp(I));
        end
    end
end
end

function [plt] = updatePlot(plt,part,stat,sim,currentStep,createPlot)
if nargin == 5
    createPlot = sim.movie;
end

for n = 1:sum(size(part.addX,2))
    line(plt.fig1.ax1,part.addX(:,n)*sim.dist.mult,part.addY(:,n)*sim.dist.mult,'Color',part.addC(n));
end

if ~createPlot
    return;
end

if sim.movie
    set(plt.fig1.plotPos, 'XData', part.x.*sim.dist.mult, 'YData', part.y.*sim.dist.mult);
    set(plt.fig1.plotVel,'XData',part.x.*sim.dist.mult,'YData',part.y.*sim.dist.mult,'UData',part.vx.*sim.dist.mult,'VData',part.vy.*sim.dist.mult);
end

title(plt.fig2.ax1,sprintf('System Temperature (Average = %.0f K)',stat.avtemp(currentStep)));
set(plt.fig2.plotTemp, 'XData', stat.time(1:currentStep)*sim.time.mult, 'YData', stat.temp(1:currentStep));
set(plt.fig2.plotAvTemp, 'XData', stat.time(1:currentStep)*sim.time.mult, 'YData', stat.avtemp(1:currentStep));

title(plt.fig4.ax1,sprintf('Current (Average = %.0f %s)',stat.avDriftCurr(currentStep)*sim.curr.mult,sim.curr.unit));
set(plt.fig4.plotCurr, 'XData', stat.time(1:currentStep)*sim.time.mult, 'YData', stat.driftCurr(1:currentStep)*sim.curr.mult);
set(plt.fig4.plotAv, 'XData', stat.time(1:currentStep)*sim.time.mult, 'YData', stat.avDriftCurr(1:currentStep)*sim.curr.mult);

delete(plt.fig3.vHist)
plt.fig3.vHist = histogram(plt.fig3.ax1,part.v*sim.speed.mult,'FaceColor','g');
plt.fig3.axis = updateAxis(plt.fig3.axis,axis(plt.fig3.ax1));
axis(plt.fig3.ax1,plt.fig3.axis);
set(plt.fig3.Vav,'Value',stat.avV(currentStep)*sim.speed.mult);
set(plt.fig3.Vrms,'Value',stat.avVrms(currentStep)*sim.speed.mult);
title(plt.fig3.ax1,sprintf('Particle Velocity Distribution at t = %f %s (Mean = %.0f %s)',stat.time(currentStep)*sim.time.mult,sim.time.unit,stat.avVrms(currentStep)*sim.speed.mult,sim.speed.unit));
end

function [out] = updateAxis(oldAxis,newAxis)
out = zeros(1,4);
out(1) = min([oldAxis(1),newAxis(1)]);
out(2) = max([oldAxis(2),newAxis(2)]);
out(3) = min([oldAxis(3),newAxis(3)]);
out(4) = max([oldAxis(4),newAxis(4)]);
end

function [] = outputCalcResults(const,sim,stat,board)
fprintf('vth = %.0f km/s\n',const.vth/1000);
fprintf('vav = %.0f km/s\n',const.vav/1000);
fprintf('Mean Free Path = %.3f nm\n',const.MFP*1e9);
fprintf('Total Sim Time = %.0f ps\n',sim.dt*sim.steps*1e12);
fprintf('Probability of Scatter = %e \n',sim.Pscat);
fprintf('Probability of Scatter = %e \n',sim.Pscat);
fprintf('Effective Mass of Electron = %e kg\n',const.meff);
fprintf('Mean Free Path (Experimental) = %f nm\n',getPlateauAverage(stat.MFP*1e9));
fprintf('Mean Time Between Collisions (Experimental) = %f ps\n',getPlateauAverage(stat.tmn*1e12));
fprintf('Mean Free Path = %.3f nm\n',mean(stat.avV)*const.tmn*1e9);
fprintf('Electric Field = <%.3e,%.3e> V/m\n',board.Ex,board.Ey);
fprintf('Electric Force = <%.3e,%.3e> N\n',board.Fx,board.Fy);
fprintf('Electric Acceleration = <%.3e,%.3e> m/s^2\n',board.ax,board.ay);
end

function [out] = getPlateauAverage(y)
y = flip(y);
s = 0;
for n = 1:length(y)
    s = s + y(n);
    if abs(s/n - y(1)) > y(1)/100
        out = s/n;
        fprintf('Average Over %d Points\n',n);
        return;
    end
end
end

function [] = plotHeatMaps(board,sim,stat)
plotParticleHeatMap(board,sim,stat);
plotTempHeatMap(board,sim,stat);
end

function [fig] = plotParticleHeatMap(board,sim,stat)
fig.h = figure;
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
axis(fig.ax1,[board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
pcolor(fig.ax1,stat.heatMap.X,stat.heatMap.Y,(stat.avHeatMap')/sim.steps);
colorbar(fig.ax1);
xlabel(fig.ax1,sprintf('x (%s)',sim.dist.unit));
ylabel(fig.ax1,sprintf('y (%s)',sim.dist.unit));
title(fig.ax1,'Particle Distribution HeatMap');
grid(fig.ax1,'on');
hold(fig.ax1, 'on')
end

function [fig] = plotTempHeatMap(board,sim,stat)
fig.h = figure;
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
axis(fig.ax1,[board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
pcolor(fig.ax1,stat.heatMap.X,stat.heatMap.Y,(stat.avTempMap')/sim.steps);
colorbar(fig.ax1);
xlabel(fig.ax1,sprintf('x (%s)',sim.dist.unit));
ylabel(fig.ax1,sprintf('y (%s)',sim.dist.unit));
title(fig.ax1,'Temperature Distribution HeatMap');
grid(fig.ax1,'on');
hold(fig.ax1, 'on')
end
