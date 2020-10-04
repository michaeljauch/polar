% This script is nearly identical to code provided by Simon Byrne and Mark Girolami, 
% the authors of 'Geodesic Monte Carlo' cited in the paper. The GMC algorithm is 
% described in full generality in Algorithm 1 of that paper. Details specific to 
% the Stiefel manifold are given in Section 4.3, while details specific to the 
% eigenmodel application are given in Section 5.3. 
 
% Read in and set up data

Y = dlmread('hoff.dat');

d = size(Y,1);
m = d; % Number of proteins in the network 
p = 3; % Rank of latent matrix of probabilities

rng(20); % Set seed

% Set step sizes for GMC algorithm
eps_u = 0.005;
eps_lambda = 0.1;
eps_c = 0.001;

% Choose the number of burn-in iterations, total iterations, and the number of 
% integration steps in the GMC algorithm. 
burn = 5000; 
N = 5000; 
T = 20;

subdiag = logical(tril(ones(d,d),-1));
Ys = Y*2-1;
Ys(1:d+1:d*d)=0;

% Set initial values 
U = randn(d,p);
[U,~] = qr(U,0);

Lambda = randn(p,1);
c = -2;

% Initialize arrays in which information from the Markov chain will be recorded. 
Us = zeros(N,d,p);
Lambdas = zeros(N,p);
cs = zeros(N,1);
lds = zeros(N,1);
accepts = zeros(N,1);

H = zeros(T,1);

% Begin GMC iterations 
for n = 1:(N + burn)
    
    % Do not save draws if we are still in the burn-in period
    if n == burn + 1
        tic
    end

    % Begin Algorithm 1 

    U_V = randn(d,p);  
    A = U'*U_V;
    U_V = U_V -0.5*U*(A+A');
    U_o = U;
    
    Lambda_v = randn(p,1);
    Lambda_o = Lambda;
    
    c_v = randn(1,1);
    c_o = c;
    
    Nu = Ys .* (U*diag(Lambda)*U' + c);
    
    ld_o = sum(log(normcdf(Nu(subdiag)))) ...
            - 0.5*sum(Lambda.^2)/m - 0.5*c^2/100;
        
    ldv_o = -0.5*sum(reshape(U_V,[],1).^2) ...
        -0.5*sum(Lambda_v.^2) -0.5*c_v^2;
        
    Nu_grad = Ys .* normpdf(Nu) ./normcdf(Nu); 
    
    U_grad = Nu_grad*bsxfun(@times,U,Lambda');
    Lambda_grad = 0.5*diag(U'*Nu_grad*U) - Lambda/m;
    c_grad = sum(Nu_grad(subdiag)) - c/100;
    
    % Initial half step 

    U_V = U_V + 0.5*eps_u*U_grad;
    A = U'*U_V;
    U_V = U_V -0.5*U*(A+A');
    
    Lambda_v = Lambda_v + 0.5*eps_lambda*Lambda_grad;
    c_v = c_v + 0.5*eps_c*c_grad;
    
    % Perform the integration steps 
    for t = 1:T
        
        
        A = U'*U_V;
        S = U_V'*U_V;
        
        tmpA = expm(-eps_u*A);
        tmpB = [U,U_V]*expm(eps_u*[A,-S;eye(p),A])*blkdiag(tmpA,tmpA);
        U = tmpB(:,1:p);
        U_V = tmpB(:,(p+1):(2*p));
        
        % re-orthonormalize if necessary
        orthcheck = norm(U'*U - eye(p),'fro');
        if orthcheck > 1e-13
            warning('re-orthonormalizing U: %s', orthcheck);
            [U,R] = qr(U,0);
            U = U*diag(sign(diag(R)));
        end
        
        Lambda = Lambda + eps_lambda*Lambda_v;
        c = c + eps_c*c_v;
        
        Nu = Ys .* (U*diag(Lambda)*U' + c);
        Nu_grad = Ys .* normpdf(Nu) ./normcdf(Nu);

        U_grad = Nu_grad*bsxfun(@times,U,Lambda');
        Lambda_grad = 0.5*diag(U'*Nu_grad*U) - Lambda/m;
        c_grad = sum(Nu_grad(subdiag)) - c/100;
    
    
        if t < T
            U_V = U_V + eps_u*U_grad;
            A = U'*U_V; 
            U_V = U_V -0.5*U*(A+A');
    
            Lambda_v = Lambda_v + eps_lambda*Lambda_grad;
            c_v = c_v + eps_c*c_grad;
        end
    end

    % Half step at the end 
    
    U_V = U_V + 0.5*eps_u*U_grad;
    A = U'*U_V; 
    U_V = U_V -0.5*U*(A+A');
    
    Lambda_v = Lambda_v + 0.5*eps_lambda*Lambda_grad;
    c_v = c_v + 0.5*eps_c*c_grad;
    
    Nu = Ys .* (U*diag(Lambda)*U' + c);
    ld = sum(log(normcdf(Nu(subdiag)))) ...
            - 0.5*sum(Lambda.^2)/m - 0.5*c^2/100;
        
    ldv = -0.5*sum(reshape(U_V,[],1).^2) ...
        -0.5*sum(Lambda_v.^2) -0.5*c_v^2;
    
    % Accept - reject step 

    if log(rand) > ld + ldv - ld_o - ldv_o
        U = U_o;
        Lambda = Lambda_o;
        c = c_o;
        if n > burn
            lds(n-burn) = ld_o;
        end
    else
        if n > burn
            accepts(n-burn) = 1;
            lds(n-burn) = ld;
        end
    end

    if n > burn
        Us(n-burn,:,:) = U;
        Lambdas(n-burn,:) = Lambda;
        cs(n-burn) = c;
    end
    
end
toc
time_after_burn = toc; 

% Save the results 
csvwrite('~/geodesic_network_eigenmodel_draws.csv', Lambdas)
csvwrite('~/geodesic_network_eigenmodel_time.csv', time_after_burn)

%mean(accepts)
%dlmwrite('~/projects/spherical/tex/network_b.dat',[(1:N)',Lambdas],' ');
plot(Lambdas(:,1))
hold on
plot(Lambdas(:,2))
hold on 
plot(Lambdas(:,3))

%plot(H);
    
