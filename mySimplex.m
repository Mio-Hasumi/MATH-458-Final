function [x_opt, f_opt] = mySimplex(c, Aeq, beq)
%   c    = [1;2;3];
%   Aeq  = [1 1 1];
%   beq  = 1;
%   [x,f] = mySimplex(c,Aeq,beq)
%   returns x = [1;0;0], f = 1

  tol = 1e-8;
  [m, n] = size(Aeq);

  % Phase I
  T = [ Aeq, eye(m), beq ];
  costRow = [ zeros(1,n), ones(1,m), 0 ];    % 0⋯0 | 1⋯1 | 0

  for i = 1:m
  costRow = costRow - T(i,:);
  end

  T = [ T; costRow ];
basis = n+1 : n+m;

  % simplex
  [T, basis] = simplexIter(T, basis, tol);

  % if optimal Phase‑I cost > 0 → infeasible
  if abs(T(end,end)) > tol
    error('LP is infeasible');
  end

  % drop artificial columns
  T(:, n+1 : n+m) = [];

  % Phase II
  costRow = [c(:)' 0];
  for i = 1:m
    j = basis(i);
    costRow = costRow - costRow(j) * T(i,:);
  end
  T(end,:) = costRow;

  % simplex
  [T, basis] = simplexIter(T, basis, tol);
  x_opt = zeros(n,1);
  for i=1:m
    if basis(i) <= n
      x_opt(basis(i)) = T(i,end);
    end
  end
  f_opt = c(:)' * x_opt;

end

%------------------------------------------------------------------------
function [T, basis] = simplexIter(T, basis, tol) % this is a helper
% Perform primal simplex
  [m1, n1] = size(T);
  while true
    cost = T(end,1:n1-1);
    [minCost, j] = min(cost);
    if minCost >= -tol
      return
    end
    col = T(1:m1-1, j);
    if all(col <= tol)
      error('LP is unbounded');
    end
    ratios = T(1:m1-1,end) ./ col;
    ratios(col <= tol) = inf;
    [~, i] = min(ratios);
    T(i,:) = T(i,:) / T(i,j);
    for k = 1:m1
      if k~=i
        T(k,:) = T(k,:) - T(k,j)*T(i,:);
      end
    end
    basis(i) = j;
  end
end
