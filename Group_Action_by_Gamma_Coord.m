function fn = Group_Action_by_Gamma_Coord(f,gamma)

[n,T] = size(f);

% fn = spline(linspace(0,1,T) , f ,gamma);
for j=1:n
      fn(j,:) = interp1(linspace(0,1,T) , f(j,:) ,gamma,'linear');
% fn(j,:) = interp1(linspace(0,1,T) , f(j,:), gamma,'linear');
end