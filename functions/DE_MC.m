function solset = DE_MC(objfun,d,lb,ub,maxGen,Npop,x2,reset)
%DE_MC Differential Evolution Markov Chain (ter Braak, 2006)
%   solset = DE_MC(objfun,d,lb,ub,maxGen,Npop,x2,reset)
%
%     objfun : A function that takes as input a vector x (which contains the
%     variables to be optimized) and the additional parameter x2. x2 serves as a
%     container for all the parameters required to run the objective function
%     that are not optimized. 
%     d : The dimensionality of the problem 
%     lb : A vector (1xd) of the lower bounds for each parameter
%     ub : A vector (1xd) of the upper bounds for each parameter
%     maxGen : The number of generations to run the optimization process
%     Npop : The number of chains to run
%     x2 : A container for all the non-optimizing parameters required to 
%     run the objective function. 
%     reset : A flag for handling boundary conditions
%             1 - Reset to boundary
%             2 - Reflect from boundry
%             3 - Periodic boundary 
% 

% Initializations
complementer = ones(d,1);
b = 1e-6; % for symmetric perturbation
% Variable declaration and initializations
X =  bsxfun(@plus,bsxfun(@times,rand(Npop,d),(ub-lb)),lb);Z = nan(Npop,d);
F = nan(Npop,1);G = nan(Npop,1);
% Solution vectors for each member and each generation
Xlist = nan(maxGen,Npop,d); 
% Objective values for each member and each generation
Flist = nan(maxGen,Npop); 
% Objective function not vectorised, use loop to set parent fitnesses.
for k=1:Npop
    F(k) = objfun(X(k,:),x2); 
end
for k=1:maxGen
    % Every 10th generation set gamma = 1
    if mod(k,10)==0
        disp(strcat('Generation : ',num2str(k)))
        gamma = 1;
    else
        gamma = 2.4/sqrt(2*d);
    end
    % Generate offspring and compute the objective function values
    for j=1:Npop
        inds = randsample([1:j-1,j+1:Npop],2);
        temp = X(j,:) + gamma*(X(inds(1),:)-X(inds(2),:))+randn*b;
        Z(j,:) = temp;
        % Flags to check boundaries 
        lflag = lb>Z(j,:); uflag = ub<Z(j,:);
        if reset==1
           % Set to bounds
           Z(j,:) = Z(j,:).*(complementer'-(lflag+uflag)) ...
               + lb.*lflag + ub.*uflag;
        elseif reset==2
           % Reflection
           Z(j,:) = Z(j,:).*(complementer'-(lflag+uflag)) + ...
           (2*lb-Z(j,:)).*lflag + (2*ub-Z(j,:)).*uflag;
           % When reflection still gives out of bounds, random initalize
           lflag = lb>Z(j,:); uflag = ub<Z(j,:);
           if sum(lflag + uflag) ~=0
               Z(j,:) = Z(j,:).*(complementer'-(lflag+uflag)) +...
               (lb + rand(1,d).*(ub-lb)).*(lflag+uflag);
           end
        elseif reset==3
            % Periodic boundary condition
            Z(j,:) = Z(j,:).*(complementer'-(lflag+uflag)) + ...
           (Z(j,:)+ub-lb).*lflag + (Z(j,:)-ub+lb).*uflag;
        end
%         G(j) = objfun(Z(j,:),x2);
    end
    % Evaluate objective functions in parallel
%     parfor j=1:Npop
    for j=1:Npop
        G(j) = objfun(Z(j,:),x2);
    end
    for j=1:Npop
%         alpha = min(1,G(j)/F(j));
        alpha = min(1,exp( (-G(j)+F(j))/2 ));
        % Solution acceptance criterion
        if rand<alpha
            X(j,:) = Z(j,:);
            F(j) = G(j);
        end
    end
    Flist(k,:) = F;
    Xlist(k,:,:) = X;
end
solset.X = X;
solset.F = F;
solset.Xlist = Xlist;
solset.Flist = Flist;