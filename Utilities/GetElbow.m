function pt = GetElbow(x)

% get vector between first and last point - this is the line
lineVec = x(end,:) - x(1);

% normalize the line vector
lineVecN = lineVec / sqrt(sum(lineVec.^2));

% find the distance from each point to the line:
% vector between all points and first point
vecFromFirst = bsxfun(@minus, x, x(1));

% Calculate the distance to the line
scalarProduct = dot(vecFromFirst, repmat(lineVecN,numel(x),1), 2);
vecFromFirstParallel = scalarProduct * lineVecN;
vecToLine = vecFromFirst - vecFromFirstParallel;

%# distance to line is the norm of vecToLine
distToLine = sqrt(sum(vecToLine.^2,2));

%# plot the distance to the line
figure('Name','distance from curve to line'), plot(distToLine)

%# now all you need is to find the maximum
[maxDist,idxOfBestPoint] = max(distToLine);

%# plot
figure, plot(curve)
hold on
plot(allCoord(idxOfBestPoint,1), allCoord(idxOfBestPoint,2), 'or')