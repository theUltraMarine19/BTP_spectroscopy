function [X_est, label, model, phi_cell] = EMfly(X, init, patch_indices, clean_indices, noise_var)
%% init
fprintf('EM for Gaussian mixture: running ... \n');
tol = 1e-3;
maxiter = 500;
J = -inf(1,maxiter);

%% Generate the phi matrices here
phi_cell = cell(numel(patch_indices), 1);
for i = 1:numel(patch_indices)
    temp_matrix = zeros(numel(clean_indices{i}), size(X,1));
    for j = 1:numel(clean_indices{i})
        temp_matrix(j, clean_indices{i}(j)) = 1;
    end
    phi_cell{i} = temp_matrix;
end

[R, X_est] = initialization(X,init, phi_cell, noise_var);


%% Generate X estimate
% X_est = rand(size(X))*255;
% scaleX = max(X(:));
% X = X./scaleX;
% X_est = X_est*scaleX;
% init_norms = sqrt(sum(X_est.^2,1));
% X_est = bsxfun(@rdivide,X_est, init_norms);
% X_est = X_est.*scaleX;

% init_norms = sqrt(sum(X.^2,1));
% X = bsxfun(@rdivide,X, init_norms);
% X = X.*scaleX;

%% main
for iter = 2:maxiter
    fprintf('iter: %d\n', iter);
    [~,label(1,:)] = max(R,[],2);
    R = R(:,unique(label));   % remove empty clusters
    model = maximization(X,X_est,R);
    [R, X_est, J(iter)] = expectation1(X,model,phi_cell, noise_var, iter);
    fprintf('Error %f\n', J(iter));
    if iter > 2
        fprintf('Error percentage : %f\n', abs(J(iter)-J(iter-1))/abs(J(iter)));
    end
    if abs(J(iter) - J(iter-1)) < tol*abs(J(iter)) 
        break; 
    end
end
end

function [R, X_est] = initialization(X, init, phi_cell, noise_var)
n = size(X,2);
if isstruct(init)  % init with a model
    [R, X_est, ~]  = expectation1(X, init, phi_cell, noise_var, 0);
elseif numel(init) == 1  % random number of components
    k = init;
    label = ceil(k*rand(1,n));
    R = full(sparse(1:n,label,1,n,k,n));
    X_est = rand(size(X))*255;
elseif all(size(init)==[1,n])  % init with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
    X_est = rand(size(X))*255;
else
    error('ERROR: init is not valid.');
end
end

function [R, X_est, J]  = expectation(X, model, phi_cell, noise_var, iter)
mu = model.mu;
Sigma = model.Sigma;
% w = model.w;

n = size(X,2);
k = size(mu,2); % num of components
R = zeros(n,k);
R1 = zeros(n,1);
X_est = zeros(size(X));
J = 0;

parfor i1 = 1:n
    j_vals = zeros(k,1);
    estimate = zeros(size(X,1),k);
    phi = phi_cell{i1};
    patch_vec = phi*X(:,i1);
    for i=1:k
        if isinf(log(det(Sigma(:,:,i))))
            j_vals(i) = inf;
            continue;
        end
        estimate(:,i) = Sigma(:,:,i)*phi'*((phi*Sigma(:,:,i)*phi' + noise_var*eye(size(phi,1)))\(patch_vec - phi*mu(:,i))) + mu(:,i);
        j_vals(i) = (patch_vec - phi*estimate(:,i))'*((noise_var*eye(size(phi,1)))\(patch_vec - phi*estimate(:,i))) + (estimate(:,i) - mu(:,i))'*(Sigma(:,:,i)\(estimate(:,i) - mu(:,i))) + log(abs(det(Sigma(:,:,i))));
    end
    [~,j_est] = min(j_vals);
    X_est(:,i1) = estimate(:,j_est);
    R1(i1,:) = j_est;
%     if mod(i1,1000) == 1
%         fprintf('patch : %d\n', i1); 
%     end
    if isinf(j_vals(j_est))
%       fprintf('jval is now infinity \n');
%     fprintf('jval : %f\n', j_vals(j_est));
    end
    J = J + j_vals(j_est);
end

for i1=1:n
    R(i1, R1(i1)) = 1;
end

% init_norms = sqrt(sum(X_est.^2,1));
% X_est = bsxfun(@rdivide,X_est, init_norms);
% X_est = X_est.*scaleX;
end

function [R, X_est, J]  = expectation1(X, model, phi_cell, noise_var, iter)
mu = model.mu;
Sigma = model.Sigma;
w = model.w;

n = size(X,2);
k = size(mu,2); % num of components
R = zeros(n,k);
X_est = zeros(size(X));
J = 0;

parfor i1 = 1:n
    j_vals = zeros(k,1);
    estimate = zeros(size(X,1),k);
    phi = phi_cell{i1};
    patch_vec = phi*X(:,i1);
    for i=1:k
        if isinf(log(det(Sigma(:,:,i))))
            j_vals(i) = inf;
            continue;
        end
        estimate(:,i) = Sigma(:,:,i)*phi'*((phi*Sigma(:,:,i)*phi' + noise_var*eye(size(phi,1)))\(patch_vec - phi*mu(:,i))) + mu(:,i);
        j_vals(i) = (patch_vec - phi*estimate(:,i))'*((noise_var*eye(size(phi,1)))\(patch_vec - phi*estimate(:,i))) + (estimate(:,i) - mu(:,i))'*(Sigma(:,:,i)\(estimate(:,i) - mu(:,i))) + log(abs(det(Sigma(:,:,i))));
    end
    [~,j_est] = min(j_vals);
    X_est(:,i1) = estimate(:,j_est);
%       if mod(i1,1000) == 1
%         fprintf('patch : %d\n', i1); 
%     end
    if isinf(j_vals(j_est))
%       fprintf('jval is now infinity \n');
%     fprintf('jval : %f\n', j_vals(j_est));
    end
    J = J + j_vals(j_est);
end

parfor j=1:k
    R(:,j) = loggausspdf(X,mu(:,j),Sigma(:,:,j));
end
R = bsxfun(@plus,R,log(w));
T = logsumexp(R,2);
R = exp(bsxfun(@minus,R,T));

end

function model = maximization(X, X_est, R)
[d,n] = size(X);
k = size(R,2); % num of components
% fprintf('k : %d\n', k);

nk = sum(R,1); % number of signals assigned to each component
w = nk/n;
mu = bsxfun(@times, X_est*R, 1./nk);

Sigma = zeros(d,d,k);
r = sqrt(R);
for i = 1:k
    Xo = bsxfun(@minus,X_est,mu(:,i)); % X - mu
    Xo = bsxfun(@times,Xo,r(:,i)');
    
%     fprintf("Num zero in XoXo' : %d\n", numel(find(Xo*Xo'==0)));
    
    Sigma(:,:,i) = Xo*Xo'/nk(i)+eye(d)*(1e-6); % to prevent non-invertibility
    if det(Sigma(:,:,i)) == 0
        fprintf('Det val : %f\n', log(abs(det(Sigma(:,:,i)))));
        fprintf('Danger!! %d\n', i);
    end
end

model.mu = mu;
model.Sigma = Sigma;
model.w = w;
end

function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;
end