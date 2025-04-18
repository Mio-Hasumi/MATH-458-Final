function [x, obj] = interiorpoint(c, Aeq, beq)
% using a primal log‑barrier + Newton’s method.
%   c   = [1;2;3];
%   Aeq = [1 1 1];
%   beq = 1;
%   [x,obj] = interiorpoint(c,Aeq,beq)

  tol       = 1e-8;     %%tolerance
  tau       = 10;       %decrease factor
  newtMax   = 50;       % max Newton steps per parameter
  [m, n]    = size(Aeq);

  x = (beq/n) * ones(n,1);
  mu = 1;

  while mu > tol
    for k = 1:newtMax
      grad = c - mu ./ x;              
      H    = diag(mu ./ (x.^2));   
      % Build
      K    = [H,        Aeq';
              Aeq, zeros(m)];
      rhs  = -[ grad; Aeq*x - beq ];
      sol  = K \ rhs;
      dx   = sol(1:n);
      dnu  = sol(n+1:end);

      % stopping criterion
      if norm(grad + Aeq'*dnu) < tol && norm(Aeq*x - beq) < tol
        break
      end

      % backtracking line‐search
      alpha = 1;
      negIdx = dx < 0;
      if any(negIdx)
        alpha = min(1, 0.99*min( -x(negIdx) ./ dx(negIdx) ));
      end
      x = x + alpha * dx;
    end
    mu = mu / tau;
  end

  obj = c'*x;
end
