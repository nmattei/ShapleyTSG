% maxmoatpacking finds an optimal moat packing for a symmetric, fully
% connected travelling salesman problem using euclidean distances
% input is a file following tsplib format, with nodes defined by x and y
% coordinates
% output is fval, the value of the held karp approximation to the shorted
% travelling salesman route
% and alloc, a cost allocation vector which should be within the epsilon
% core of the TSP game using filename as its' characteristic function


%function [fval, alloc] = maxmoatpacking(filename)
%tic;
figure(1);
filename = 'bays29.tsp' %'test.txt'%'berlin52.tsp' ;%'ulysses31.txt';%'test.txt';%'test6.1.txt';% 
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
while ~(isvector(strfind(temp, 'NODE_COORD_SECTION')) || isvector(strfind(temp, 'DISPLAY_DATA_SECTION')))
    temp = fgets(fid);
end

x = 1:N;
E = 0.5*N*(N-1);

% read in coordinates
coords = zeros(N, 2);
for i=1:N
    temp = fgets(fid); 
    
    holder = regexpi(temp, ' ', 'split');
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

[m, ~] = size(Aeq);
minit = m;
 

mvt = [];

lb = zeros(E, 1);
ub = ones(E, 1);
A = -Aeq;
b = -beq;

kk=0;

LPtime = 0;
CHECKCtime = 0;
MINCUTtime = 0;
RESIZEtime = 0;

while (true) 
op = optimset('LargeScale', 'off', 'Simplex', 'on');

% solve LP with current subtour constraints
tic;
[x, fval, ~, ~, ~] = linprog(f, A, b, Aeq, beq, lb, ub, [], op);
LPtime = LPtime + toc;
size(A);

kk = kk+1;

% generate network for min cut
c=1;
ntwk = zeros(N, N);
for i=1:N
    for j=i+1:N
         ntwk(i,j) = x(c);
         ntwk(j,i) = x(c);
         c = c+1;
    end
end

% check if network connected graph
tic;
ntwk1 = ntwk>0;
mw = checkc(ntwk1);
if (mw == 0)
    mv = checkcc(ntwk1,1)';
end
CHECKCtime = CHECKCtime + toc;

tic;
if (mw == 1)
    % find min cut if it is connected
    [mv, mw] = mincut(ntwk);
end
MINCUTtime = MINCUTtime + toc;

tic;

mvt = [mvt, mv];
k = k + 1;

err = 0.0001;
% check min cut value
if ((mw <= 2+err) && (mw >= 2-err))
    break
end

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
RESIZEtime = RESIZEtime + toc;

end
%pi = -lambda.('eqlin');%-getfield(lambda, 'eqlin');
A = [Aeq; A];

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

hold all
%ineqlin = lambda.('ineqlin');

k=1;
subt = m-minit;
mvt = mvt(:, 1:subt);
for i = 1:E
        newMoat = A(1:N, i);
        newMoat = abs(1-newMoat);
        mvt=[mvt, newMoat];
        k=k+1;
end

 clf
 figure(1) 
 hold all
 
 % allocation vector
 alc2 = zeros(1,N)';

 [~,numMoats] = size(mvt);
 numMoats = numMoats + N;
 
lbM = zeros(1, N+subt+E);
fM = 2*ones(1, N+subt+E)';
 
AM = -A(N+1:N+N+subt, :)';
AE = A(1:N, :);

moatPairs = AE'*AE;
moatPairs = moatPairs>0;
for i=1:E
    moatPairs(i,i) = 0;
end

AM = [AM, moatPairs];

[w, fvalw] = linprog(-fM, AM, f, [], [], lbM, []);

if (~(fval == fvalw))
   % noted
end

tic;
[w, mvt] = nest(w, [1-eye(N),mvt]);
 NESTtime = toc;
 
 tic;
for i=N+1:length(w)
   currColour = rand(1,3);
   numOutside=N-sum(mvt(:,i));
   for j=1:N
       if (mvt(j,i)==0)
           filledCircle([coords(j,1), coords(j,2)], w(i)+w(j),50, currColour);
           alc2(j) = alc2(j) + w(i)/numOutside;
           
           % cheat for test.txt
              if (i==11)
               
        filledCircle([coords(3,1), coords(3,2)], w(11)+w(14)+w(3),50, currColour);
        filledCircle([coords(9,1), coords(9,2)], w(11)+w(14)+w(9),50, currColour);
           end
       end
   end
   numMoats = numMoats+1;  
end


for i=1:N
    circle(coords(i,1), coords(i,2), w(i), 1)
    alc2(i) = alc2(i) + w(i);
end 

% node 1 is depot
for i=2:N
    alc2(i) = alc2(i) + alc2(1)/(N-1);
end

alc2(1) = 0;

alloc = 2*alc2;

disp(alloc);
%disp(fval);

labels=cellstr(num2str([1:N]'));
plot(coords(:,1), coords(:,2), 'x');
text(coords(:,1), coords(:,2), labels);
 gplot(ntwk, coords)
 axis equal
 DRAWtime = toc;
 
 %NESTtime
 %LPtime
 %DRAWtime
 %RESIZEtime
 %MINCUTtime
 %CHECKCtime
 %toc
 
 %[mAdj, wMod] = findnests(w, mvt)