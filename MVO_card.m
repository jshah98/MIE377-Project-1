function  x = MVO_card(mu, Q, targetRet, card)

    % Use this function to construct your cardinality-constrainted MVO
    % portfolio subject to the target return, with short-selling
    % disallowed.
    %
    % You may use Gurobi or any other optimizer that is able to solve mixed
    % integer programs (MIPs). Just be sure to comment on your code to
    % (briefly) explain your procedure.

    % Find the total number of assets
    n = size(Q,1);

    B = [-eye(n) lBuy*eye(n); eye(n) -uBuy*eye(n)];
    A = [B;-1*transpose(mu) ];
    b = [zeros(2*n,1);-targetRet];
    Aeq = [ones(1,n) zeros(1,n); zeros(1,n) ones(1,n)];
    beq = [1;card];
    varTypes = [repmat('C', n, 1);repmat('B', n, 1)];
    lb = zeros(2*n,1);
    ub = ones(2*n,1);

    clear model;

    model.Q = sparse(2*Q);
    model.obj = zeros(1,2*n);
    model.A = [sparse(A);sparse(Aeq)];
    model.rhs = full([b; beq]);
    model.sense = [ repmat('<', (2*n + 1), 1) ; repmat('=', 2, 1) ];
    model.vtype = varTypes;
    model.lb = lb;
    model.ub = ub;

    clear params;
    params.TimeLimit = 100;
    params.OutputFlag = 0;

    x = gurobi(model,params);

end
