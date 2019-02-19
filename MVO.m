function  x = MVO(mu, Q, targetRet, varargin)

    % Use this function to construct your MVO portfolio subject to the
    % target return, with short-selling disallowed.
    %
    % You may use quadprog, Gurobi, or any other optimizer you are familiar
    % with. Just be sure to comment on your code to (briefly) explain your
    % procedure.

    % Find the total number of assets
    n = size(Q,1);
    f = zeros(n, 1);
    
    % Increasing the tolerance of 'quadprog'
    options = optimoptions( 'quadprog', 'TolFun', 1e-9 );
    % quadprog (below) minimizes 0.5*x' * Q *x s.t. -1*mu*x <= targetRet, [col of ones] * x = 1, and 0 <= xi <= 1

    x = quadprog(Q, f, -1*mu, targetRet, ones(1, n), 1, zeros(1,n), ones(1,n),[],options)
end
