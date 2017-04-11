function h = distance_(x,y, dir)
%DISTANCE_ Determine distances between locations
%
%   This function produces a matrix that describes the
%   distances between two sets of locations.
%
%INPUT PARAMETERS
%   x   - location coordinates for data set #1 [n1 x D]
%   y   - location coordinates for data set #2 [n2 x D]
%   dir - if given, specifies the component direction to use 
%OUTPUT PARAMETERS
%   h - distance between points in x from points in y [n1 x n2]
%
%EXAMPLE:
%   x=[0,0; 5,0; 5,5]
%   y=[1,1; 5,5]
%   h = distance_(x,y)
%RESULT:
%    1.4142    7.0711
%    4.1231    5.0000
%    5.6569         0

%EXTERNAL FUNCTIONS CALLED: none
%REVISION HISTORY:
%   pkk, 5/13/97
%   tae, 6/26/98 (minor changes)
%   jlm, 4/11/2017 (add component direction)

if nargin < 3, dir = []; end

[n1,D] = size(x);
[n2,D2] = size(y);

if D~=D2
   error('ERROR in DISTANCE_: locations must have same number of dimensions (columns)')
end
h = zeros(n1,n2);

if isempty(dir)
    for id = 1:D
       h = h + (x(:,id)*ones(1, n2)-ones(n1,1)*y(:,id)').^2;
    end
else
    h = (x(:,dir)*ones(1, n2)-ones(n1,1)*y(:,dir)').^2;
end
h = sqrt(h);

end