function  [mu, Q] = PCA(returns, varargin)

    % Use this function to perform the principal component analysis. After
    % performing the PCA, use the principal components and the eigenvalues
    % to construct a factor model and estimate mu and Q.
    %
    % Note: PCA only requires the asset returns, it does not need the
    % factor returns.

    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------








    % mu =          % n x 1 vector of asset exp. returns
    % Q  =          % n x n asset covariance matrix
    %----------------------------------------------------------------------
    cov_mat = cov(returns); % creates a covariance matrix of returns
    [V,D] = eig(cov_mat) % V -> matrix of eigenvectors, D -> diag(eigenvalues corresponding to vectors in V)
    [m,n] = size(returns) % gets how big f should be if we don't have K
    f = ones(m, n);

     V_inv = inv(V);
    for r = 1:m
      for c = 1:n
        f(r, c) =returns(r, :)*V_inv(:, c);
      end
    end
    [mu, Q] = FF(returns, f);
end
