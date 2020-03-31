%function [fval, alloc] = maxmoatpacking(filename)
%format compact
figure(1);
filename = 'bays29.tsp'
%filename = 'ulysses22.tsp';%'test.txt';%'test6.1.txt';% 
fid = fopen(filename);
temp = fgets(fid);
while ~(isvector(strfind(temp, 'DIMENSION')))
    temp = fgets(fid);
end
clf
N = str2double(temp(11:length(temp)-1)); % for Charles format
if (~(N>0)) % for tsplib format 
    N = str2double(temp(12:length(temp)-1)); 
end

while ~(isvector(strfind(temp, 'EDGE_WEIGHT_SECTION')))
   temp = fgets(fid);    
end
dd = zeros(N, N);
for i=1:N
        temp = fgets(fid);
       % temp = strtrim(temp);
         holder1 = regexpi(temp, ' ', 'split');
        %# find empty cells
        emptyCells = cellfun(@isempty,holder1);
         %# remove empty cells
        holder1(emptyCells) = [];
        for j=1:N
            dd(i,j) = str2double(holder1(j));
        end
end


while ~(isvector(strfind(temp, 'NODE_COORD_SECTION')) || isvector(strfind(temp, 'DISPLAY_DATA_SECTION')))
    temp = fgets(fid);
end

x = 1:N;
E = 0.5*N*(N-1);

% read in coordinates
coords = zeros(N, 2);
for i=1:N
    temp = fgets(fid); 
    temp = strtrim(temp);
    holder = regexpi(temp, ' ', 'split');
  %# find empty cells

emptyCells = cellfun(@isempty,holder);
%# remove empty cells
holder(emptyCells) = [];
    x = cell2mat( holder(2));
    y = cell2mat(holder(3));
    coords(i, 1) = str2double(x);
    coords(i, 2) = str2double(y); 
    
end

% create distance matrix and vector
d = zeros(N, N);
f = zeros(1,E);

c = 1;
for i=1:N
    for j=i+1:N
           dist = ((coords(i,1)-coords(j,1))^2 + (coords(i,2)-coords(j,2))^2)^0.5;
           d(i,j) = dist;
           d(j,i) = dist;
        f(c) = d(i,j);
        c = c+1;
    end
end

% generate linear system
Aeq = zeros(N, E);

c = 1;
for i=1:N
    for j=c:c+N-i-1
            Aeq(i,j)=1;
            c = c + 1;
    end
end

c = 1;
for i=1:N-1
    for j=1:N-i
        Aeq(i+j, c+j-1)=1;
    end
    c = c+N-i;
end

beq = 2*ones(N, 1);

k = 0;

Aub = [eye(E), eye(E)];
bub = ones(E, 1);

%A = [Aeq, zeros(N, E); Aub];
%b = [beq; bub];
A = [];
b = [];
[m, n] = size(Aeq);
minit = m;
f = [f]; %zeros(1, E)];

ineqlinS = [];
mvt = [];

lb = zeros(E, 1);
ub = ones(E, 1);
A = -Aeq;
b = -beq;
kk=0
while (true) 
op = optimset('LargeScale', 'off', 'Simplex', 'on');
[x, fval, exitflag, output, lambda] = linprog(f, A, b, Aeq, beq, lb, ub, [], op)

kk = kk+1;
%[fval, x, pi, ~, exitflag] = fullrsm(A, b, f', m, n);

c=1;
ntwk = zeros(N, N);
for i=1:N
    for j=i+1:N
         ntwk(i,j) = x(c);
         ntwk(j,i) = x(c);
         c = c+1;
    end
end

% check if ntwk connected graph
mw = checkc(ntwk);
if (mw == 0)
    mv = checkcc(ntwk,1)';
end
if (mw == 1)
    % find min cut if is connected
    [mv, mw] = mincut(ntwk);
end

mvt = [mvt, mv];
k = k + 1;

err = 0.0001;
% check min cut value
if ((mw <= 2+err) && (mw >= 2-err))
   % disp('done')
    break
end
%ineqlin = getfield(lambda, 'ineqlin');
%if (length(ineqlin) > N)
 %   ineqlinS = [ineqlinS, ineqlin(N+1)];
%end
% add subtour constraint
newrow = zeros(1, E);
for i=1:N
    
    if (mv(i)==0)
        % on outside of this 'subtour moat', remove edges to other nodes in
        % set 0
        newrow = abs(newrow-Aeq(i, :));
    end
end

A = [A; -newrow];
b = [b; -2];
m = m + 1;



end
pi = -getfield(lambda, 'eqlin');
A = [Aeq; A];
b = [beq; b];

% final adjacency matrix of approximated tour - not needed
c=1;
xe = x(1:E);
ntwk = zeros(N, N);
for i=1:N
    for j=i+1:N
         ntwk(i,j) = xe(c);
         ntwk(j,i) = xe(c);
         c = c+1;
    end
end

 clf
% plots and allocations, plots not needed
hold all
colours = [0.3 0.1 0.6; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1];
numMoats = 0;

alloc = zeros(N,1);


ineqlin = getfield(lambda, 'ineqlin');
% moats from added subtour constraints


ubduals = pi(N+1:size(pi)-(m-minit)-1);
At = A';

ubduals = getfield(lambda, 'upper');
lbduals = getfield(lambda, 'lower');
offset = zeros(N,1);
offsetL = zeros(N,1);

%for i=1:N
%    rowAt = At(:,i);
%    offsetTemp = rowAt*ubduals;
%    offsetLow = rowAt*lbduals;
%    if (offsetTemp ~= 0)
%        numOutside = N-sum(At(i, 1:N));
%        offset(1) = offsetTemp;
%        if (At(i, 1) == 0) % 0 is inside if node1 is 0
%            for j=2:N
%                if (At(i,j)==1)
%                    alloc(j) = alloc(j) - offset(j)*2/(N-numOutside);
%                   offset(j) = offset(j) + offsetTemp;          
%                    offsetL(j) = offsetL(j) + offsetLow;  
%                end
%            end
 %       else % 1 is inside if node1 is 1
 %           for j=2:N
 %                if (At(i,j)==1)
 %                    alloc(j) = alloc(j) - offset(j)*2/numOutside;
 %                    offset(j) = offset(j) + offsetTemp;
 %                    offsetL(j) = offsetL(j) + offsetLow;  
 %               end
 %           end
 %       end
 %   end
moatW = ineqlin;
At = A';
offset = ubduals' * At;
offset = offset(1:N);
k=1;

mvt = mvt(:, 1:m-minit);
for i = 1:E
  %  if (ubduals(i) ~= 0)
        newMoat = A(1:N, i);
     %   if (newMoat(1) == 0)
            newMoat = abs(1-newMoat);
   %     end
        mvt=[mvt, newMoat];
        moatW(length(ineqlin)+k) = ubduals(i)*0.5;
        k=k+1;
   % end
end

[~,moatN] = size(mvt);
colours = colours(1:6, :);

for i=1:N
    if (offset(i) ~= 0)
        pi(i) = pi(i) - 0.5* offset(i);
        
        %        circle(coords(i,1), coords(i,2), pi(i)-offset(i)*0.5, 1)
    %    filledCircle([coords(i,1),coords(i,2)], pi(i), 50, 'y');
   %      circle(coords(i,1), coords(i,2), pi(i)+0.5*offset(i), 1);
    %   % filledCircle([coords(i,1),coords(i,2)], pi(i), 50, 'y');
   % %   circle(coords(i,1), coords(i,2), pi(i), 1)
         numMoats = numMoats+1;
    end
end

for i=1:moatN
   currColour = colours(mod(i,6)+1, :);
   numOutside=N-sum(mvt(:,i));
   for j=1:N
       if (mvt(j,i)==0)
           filledCircle([coords(j,1), coords(j,2)], pi(j)+moatW(i+N),50, currColour);
          alloc(j) = alloc(j) +  moatW(i+N)/numOutside;
            
       end
   end
   numMoats = numMoats+1;  
end


% moats from slack var duals


% indv moats, adjusted
tol = 0.0001;
for i=1:N
  %  if (pi(i)+0.5*offset(i) >= tol)
    circle(coords(i,1), coords(i,2), pi(i), 1)
    alloc(i) = alloc(i) + pi(i);
    numMoats = numMoats+1;
  %  end
end 

for i=2:N
    alloc(i) = alloc(i) + alloc(1)/(N-1);
end
alloc(1) = 0;
alloc = alloc*2;

moats = zeros(numMoats,1);
for i=1:N
    moats(i) = pi(i)+0.5*offset(i);
end

%[~, w]  = size(mvt);
%???????????????????????????????????????? ??????????
%for i=E+N+1:m
%   currCol = mvt(:,i-E-N);
%    numOutside=N-sum(mvt(:, i-E-N));
%   for j=1:N
%       if (currCol(j) == 0)
%           alloc(j) = alloc(j) + 2*pi(i)/numOutside;
%       end
%   end
%end

for i=2:N
 %   alloc(i) = alloc(i) + 2*moats(i);
   % alloc(i) = alloc(i) + 2*moats(1)/(N-1);
end

numMoats
subtoursElim = k-1
scaled = (alloc/sum(alloc))*fval

labels=cellstr(num2str([1:N]'));
plot(coords(:,1), coords(:,2), 'x');
text(coords(:,1), coords(:,2), labels);
 gplot(ntwk, coords)
 axis equal
 
 
 figure(2)
 
 alc2 = zeros(1,N)';
 [~, wdth] = size(mvt);
 %mvt2=mvt;
mvt2 = mvt(:, 1:(wdth-(sum(ubduals>0))));
 clf
 hold all
 [~,numMoats] = size(mvt);
 numMoats = numMoats + N;
 totM = numMoats + N + (m-minit);
 lbM = zeros(1, N+(m-minit)+E);
 fM = 2*ones(1, N+(m-minit)+E)';%+E)';
 
AM = -A(N+1:N+N+(m-minit), :)';
AE = A(1:N, :)
%AM2 = [AM, [AE; zeros(E-N,E)]]

subt = m-minit;
moatPairs = AE'*AE;
moatPairs = moatPairs>0;
for i=1:E
    moatPairs(i,i) = 0;
end

AMo = AM;
fMo = fM(1:N+subt);
AM = [AM, moatPairs];

[w, fvalw] = linprog(-fM, AM, f, [], [], lbM, []);
[ wo, fvalwo] = linprog(-fMo, AMo, f,  [], [], lbM(1:E), []);
wi = w(1:N);
wb = w(N+1:N+E+subt);
[w, mvt] = nest(w, [1-eye(N),mvt]);
[mvtH, mvtW] = size(mvt);
%w = [wi; wb2];
lw = length(w);
wOrig = w;
ws = w(1:N);
%[mAdj, mvte] = findnests(w(N+1:lw), mvt(:,N+1:lw));
%w = [ws; we];
%mvt = [1-eye(N), mvte];
currColour = rand(1,3);
 
 
 
for i=N+1:length(w)
   %currColour = colours(mod(i,6)+1, :);
   currColour = rand(1,3);
   numOutside=N-sum(mvt(:,i));
   for j=1:N
       if (mvt(j,i)==0)
           filledCircle([coords(j,1), coords(j,2)], w(i)+w(j),50, currColour);
            alc2(j) = alc2(j) + w(i)/numOutside;
          % if (i==11)
               
        %filledCircle([coords(3,1), coords(3,2)], w(11)+w(14)+w(3),50, currColour);
        %filledCircle([coords(9,1), coords(9,2)], w(11)+w(14)+w(9),50, currColour);
       %    end
       end
   end
   numMoats = numMoats+1;  
end


for i=1:N
    
    circle(coords(i,1), coords(i,2), w(i), 1)
    
    alc2(i) = alc2(i) + w(i);
    %for j=1:N
      %  if (i~=j)
     %       alc2(j) = alc2(j) + w(i)/(N-1);
    %    end
   % end
end 

for i=2:N
    alc2(i) = alc2(i) + alc2(1)/(N-1);
end
alc2(1) = 0;

alc2=2*alc2;
allocB = alloc;
alloc = alc2

%[alloc, alc2]

labels=cellstr(num2str([1:N]'));
plot(coords(:,1), coords(:,2), 'x');
text(coords(:,1), coords(:,2), labels);
 gplot(ntwk, coords)
 axis equal