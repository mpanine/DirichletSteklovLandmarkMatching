function [ Energies, dEnergies ] = compute_full( Src, Tar, energy_settings )
%COMPUTE_FULL Computes the energy term Energy for the functional map evaluation and its
%derivative dEnergy in the case of full Fmaps.
    assert(isa(energy_settings, 'energy.Settings'), ...
        'energy_settings must be of type energy.Settings or energy.PartialSettings.');
    
    Energies = {};
    dEnergies = {};


    numBasisFun_Src = energy_settings.basis_settings_Src.numBasisFun;
    numBasisFun_Tar = energy_settings.basis_settings_Tar.numBasisFun;

    basis_Src = Src.BASIS.basis(:,1:numBasisFun_Src);
    basis_Tar = Tar.BASIS.basis(:,1:numBasisFun_Tar);
    basis_inverse_Src = Src.BASIS.basis_inverse(1:numBasisFun_Src,:);
    basis_inverse_Tar = Tar.BASIS.basis_inverse(1:numBasisFun_Tar,:);
    
    a = energy_settings.weight_descriptorPreservation;
    b = energy_settings.weight_descriptorCommutativity;
    c = energy_settings.weight_LBCommutativity;
    alpha = energy_settings.weight_orientationTerm;
    
    if energy_settings.useLegacyNormalization
        [C_eye,~,~] = mat_projection(eye(numBasisFun_Tar, numBasisFun_Src));
        C_eye = reshape(C_eye,...
                        [numBasisFun_Tar*numBasisFun_Src, 1]);
%         C_eye = reshape(eye(numBasisFun_Tar, numBasisFun_Src),...
%                         [numBasisFun_Tar*numBasisFun_Src, 1]);
    end
                
    if energy_settings.computeDescriptorPreservation || ...
       energy_settings.computeDescriptorCommutativity || ...
       energy_settings.computeOrientationTerm
        assert(isfield(Src, 'DESCRIPTORS') && isfield(Tar, 'DESCRIPTORS'), ...
            'The Source and Target shapes must have descriptors.');

        fct_Src = Src.DESCRIPTORS.fct;
        fct_Tar = Tar.DESCRIPTORS.fct;
        assert(size(fct_Src,2)==size(fct_Tar,2), ...
            'The number of descriptors on Source and Target must be equal.');

        % Normalization of the descriptor functions
        no_Src = sqrt(diag(fct_Src'*Src.SHAPE.A*fct_Src))';
        fct_Src = fct_Src ./ repmat(no_Src, [Src.SHAPE.nv,1]);
        no_Tar = sqrt(diag(fct_Tar'*Tar.SHAPE.A*fct_Tar))';
        fct_Tar = fct_Tar ./ repmat(no_Tar, [Tar.SHAPE.nv,1]);
    end

    if energy_settings.computeDescriptorPreservation
        % ---
        % Descriptors preservation term
        Fct_src = basis_inverse_Src*fct_Src;
        Fct_tar = basis_inverse_Tar*fct_Tar;

        if energy_settings.useOperatorNormalization
            max_norm = max(norm(Fct_src,'fro'), norm(Fct_tar,'fro'));
            Fct_src = Fct_src/max_norm;
            Fct_tar = Fct_tar/max_norm;
        end

        E_descriptorPreservation = @(FF) sum(sum(...
            (reshape(FF, [numBasisFun_Tar,numBasisFun_Src])*Fct_src - Fct_tar).^2))/2;
        dE_descriptorPreservation = @(FF) vec(...
            (reshape(FF, [numBasisFun_Tar,numBasisFun_Src])*Fct_src - Fct_tar)*Fct_src');
        
        if energy_settings.useLegacyNormalization
            a_nor = a/E_descriptorPreservation(C_eye);
            if ~isinf(a_nor)
                a = a_nor;
            end
        end
        
        E_descriptorPreservation = @(FF) a*E_descriptorPreservation(FF);
        dE_descriptorPreservation = @(FF) a*dE_descriptorPreservation(FF);
        % ---

        Energies{end+1} = E_descriptorPreservation;
        dEnergies{end+1} = dE_descriptorPreservation;
    end

    if energy_settings.computeDescriptorCommutativity
        % ---
        % Descriptor commutativity term
        numFct = size(fct_Src,2);
        OpSrc = cell(numFct,1);
        OpTar = cell(numFct,1);
        for i = 1:numFct
            OpSrc{i} = basis_inverse_Src*(...
                repmat(fct_Src(:,i), [1,numBasisFun_Src]).*basis_Src);
            OpTar{i} = basis_inverse_Tar*(repmat(...
                fct_Tar(:,i), [1,numBasisFun_Tar]).*basis_Tar);

            if energy_settings.useOperatorNormalization
                max_norm = max(norm(OpSrc{i},'fro'), norm(OpTar{i},'fro'));
                OpSrc{i} = OpSrc{i}/max_norm;
                OpTar{i} = OpTar{i}/max_norm;
            end
        end

        E_descriptorCommutativity = @(F) sum(...
            cell2mat(cellfun(@(X,Y) sum(sum((...
            X*reshape(F, [numBasisFun_Tar,numBasisFun_Src]) - ...
            reshape(F, [numBasisFun_Tar,numBasisFun_Src])*Y).^2))...
            , OpTar', OpSrc', 'UniformOutput', false)), 2)/2;
        dE_descriptorCommutativity = @(F) sum(...
            cell2mat(cellfun(@(X,Y) vec(...
            X'*(X*reshape(F, [numBasisFun_Tar,numBasisFun_Src]) - ...
            reshape(F, [numBasisFun_Tar,numBasisFun_Src])*Y) - ...
            (X*reshape(F, [numBasisFun_Tar,numBasisFun_Src]) - ...
            reshape(F, [numBasisFun_Tar,numBasisFun_Src])*Y)*Y')...
            , OpTar', OpSrc', 'UniformOutput', false)), 2);
        % ---
        
        if energy_settings.useLegacyNormalization
            b_nor = b/E_descriptorCommutativity(C_eye);
            if ~isinf(b_nor)
                b = b_nor;
            end
        end
        
        E_descriptorCommutativity = @(F) b*E_descriptorCommutativity(F);
        dE_descriptorCommutativity = @(F) b*dE_descriptorCommutativity(F);
        
        Energies{end+1} = E_descriptorCommutativity;
        dEnergies{end+1} = dE_descriptorCommutativity;
    end

    if energy_settings.computeLBCommutativity
        if energy_settings.useResolvent
            if ~strcmp(Src.BASIS.settings.name,'MHB') ||...
               ~strcmp(Tar.BASIS.settings.name,'MHB')
                throw(MException('ENERGYCOMPUTE_FULLERROR:incorrectInputBasis',...
                    'The resolvent can only be used with a basis of type MHB defined on the Source and Target shapes.'));
            end
            
            Ev1 = Src.BASIS.eigenvalues;
            Ev2 = Tar.BASIS.eigenvalues;
            
            numEigsSrc = size(Ev1, 1);
            numEigsTar = size(Ev2, 1);

            max_ev = max(max(Ev1),max(Ev2));
            Ev1 = Ev1./max_ev;
            Ev2 = Ev2./max_ev;

            sigma = 0.5;
            Ev1 = Ev1.^sigma;
            Ev2 = Ev2.^sigma;

            % real part
            x1 = Ev2./(1 + Ev2.^2);
            x2 = Ev1'./(1+ Ev1'.^2);
            Dlb1 = (repmat(x1, [1,numEigsSrc]) - repmat(x2, [numEigsTar,1])).^2;

            % imaginary part
            x1 = 1./(1 + Ev2.^2);
            x2 = 1./(1+ Ev1'.^2);
            Dlb2 = (repmat(x1, [1,numEigsSrc]) - repmat(x2, [numEigsTar,1])).^2;

            if energy_settings.useOperatorNormalization
                max_norm = max(norm(Dlb1,'fro'), norm(Dlb2,'fro'));
                Dlb1 = Dlb1/max_norm;
                Dlb2 = Dlb2/max_norm;
            end
            
            % equal weight
            Dlb = (Dlb1 + Dlb2);
            Dlb = Dlb/norm(Dlb,'fro');

            E_LBCommutativity = @(F) sum(sum((reshape(F, [numBasisFun_Tar,numBasisFun_Src]).^2 .* Dlb)/2));
            dE_LBCommutativity = @(F) reshape(reshape(F, [numBasisFun_Tar,numBasisFun_Src]).*Dlb,[],1);

        else
%             Ev1 = Src.BASIS.eigenvalues;
%             Ev2 = Tar.BASIS.eigenvalues;
%             
%             numEigsSrc = size(Ev1, 1);
%             numEigsTar = size(Ev2, 1);
% 
%             max_ev = max(max(Ev1),max(Ev2));
%             Ev1 = Ev1./max_ev;
%             Ev2 = Ev2./max_ev;
%             
%             Dlb = (repmat(Ev2, [1,numEigsSrc]) - repmat(Ev1', [numEigsTar,1])).^2;
%             Dlb = Dlb/norm(Dlb, 'fro');
%             
%             E_LBCommutativity = @(C) sum(sum((reshape(C,[numEigsTar,numEigsSrc]).^2 .* Dlb)/2));
%             dE_LBCommutativity =@(C) reshape(reshape(C,[numEigsTar,numEigsSrc]).*Dlb,[],1);
            % Laplace-Beltrami operator computation on Source and Target Shapes
            OpSrc = cell(1,1);
            OpTar = cell(1,1);
            OpSrc{1} = ((Src.SHAPE.A*basis_Src)\(Src.SHAPE.W*basis_Src));
            OpTar{1} = ((Tar.SHAPE.A*basis_Tar)\(Tar.SHAPE.W*basis_Tar));

            if energy_settings.useLegacyNormalization
                max_norm = max(max(diag(OpSrc{1})),max(diag(OpTar{1})));
                OpSrc_nor = OpSrc{1}/max_norm;
                OpTar_nor = OpTar{1}/max_norm;
                
                if ~isinf(OpSrc_nor)
                    OpSrc{1} = OpSrc_nor;
                end
                if ~isinf(OpTar_nor)
                    OpTar{1} = OpTar_nor;
                end
            end

            % ---
            % Laplace-Beltrami operator commutativity preservation term
            E_LBCommutativity = @(F) sum(...
                cell2mat(cellfun(@(X,Y) sum(sum((...
                X*reshape(F, [numBasisFun_Tar,numBasisFun_Src]) - ...
                reshape(F, [numBasisFun_Tar,numBasisFun_Src])*Y).^2))...
                , OpTar', OpSrc', 'UniformOutput', false)), 2)/2;
            dE_LBCommutativity = @(F) sum(...
                cell2mat(cellfun(@(X,Y) vec(...
                X'*(X*reshape(F, [numBasisFun_Tar,numBasisFun_Src]) - ...
                reshape(F, [numBasisFun_Tar,numBasisFun_Src])*Y) - ...
                (X*reshape(F, [numBasisFun_Tar,numBasisFun_Src]) - ...
                reshape(F, [numBasisFun_Tar,numBasisFun_Src])*Y)*Y')...
                , OpTar', OpSrc', 'UniformOutput', false)), 2);
            % ---
        end

        
        if energy_settings.useLegacyNormalization
            if strcmp(Src.BASIS.settings.name,'MHB') && strcmp(Tar.BASIS.settings.name,'MHB')
                Ev1 = Src.BASIS.eigenvalues;
                Ev2 = Tar.BASIS.eigenvalues;
            else
                Ev1 = Src.SHAPE.laplacianBasis.eigenvalues(1:numBasisFun_Src);
                Ev2 = Tar.SHAPE.laplacianBasis.eigenvalues(1:numBasisFun_Tar);
            end;
            max_ev = max(max(Ev1),max(Ev2));
            Ev1 = Ev1./max_ev;
            Ev2 = Ev2./max_ev;
            
            if ~energy_settings.useResolvent
                Dlb = (repmat(Ev2, [1,numBasisFun_Src]) - repmat(Ev1', [numBasisFun_Tar,1])).^2;
                c = c/norm(Dlb,'fro');
            end

        
            if energy_settings.useResolvent
                sigma = 0.5;
                Ev1 = Ev1.^sigma;
                Ev2 = Ev2.^sigma;
            end
            Dlb1 = (repmat(Ev2, [1,numBasisFun_Src]) - repmat(Ev1', [numBasisFun_Tar,1])).^2;
            Dlb1 = Dlb1/norm(Dlb1, 'fro'); % make the weight consistent across different masks
        
            
            if max(abs(diag(Dlb1))) >= 0.000001*min(min(abs(Ev1(2:end))),min(abs(Ev2(2:end)))) % Workaround for isospectral case (i.e. using the same shape as source and target)
                c = c/(sum(sum((reshape(C_eye, [numBasisFun_Tar, numBasisFun_Src]).^2 .* Dlb1)/2)));
            end
        end
            
        E_LBCommutativity = @(F) c*E_LBCommutativity(F);
        dE_LBCommutativity = @(F) c*dE_LBCommutativity(F);
        
        Energies{end+1} = E_LBCommutativity;
        dEnergies{end+1} = dE_LBCommutativity;
    end

    if energy_settings.computeOrientationTerm
        assert(isfield(Src.SHAPE, 'normals_face'), 'Source shape has to have a normals_face field (set computeNormals to true).');
        assert(isfield(Tar.SHAPE, 'normals_face'), 'Target shape has to have a normals_face field (set computeNormals to true).');

        % Orientation-preserving Operators
        compute_all_OrientationOp = @(S,B,fct) ...
            cellfun(@(f) OrientationOp(S,B,f),mat2cell(fct,size(fct,1),ones(size(fct,2),1)),'un',0);
        F11_all = compute_all_OrientationOp(Src.SHAPE,Src.BASIS.basis,fct_Src);
        F22_all = compute_all_OrientationOp(Tar.SHAPE,Tar.BASIS.basis,fct_Tar);
        
        if energy_settings.useOperatorNormalization
            max_norm = num2cell(cellfun( @(F11,F22) max(norm(F11,'fro'), norm(F22,'fro')),...
                                F11_all, F22_all));
            F11_all = cellfun( @(F11, max_n) F11/max_n,...
                                F11_all, max_norm,'un',0 );
            F22_all = cellfun( @(F22, max_n) F22/max_n,...
                                F22_all, max_norm,'un',0 );
        end

        % Orientation term
        E_orientationTerm = @(F) sum(...
                        cellfun(@(F11,F22) ...
                        0.5*norm(...
                        reshape(F, [numBasisFun_Tar,numBasisFun_Src])*F11 - ...
                        F22*reshape(F, [numBasisFun_Tar,numBasisFun_Src]), 'fro')^2,...
                        F11_all, F22_all));
        dE_orientationTerm = @(F) sum(...
                        cell2mat(cellfun(@(F11,F22) ...
                        reshape(reshape(F, [numBasisFun_Tar,numBasisFun_Src])*(F11*F11') - ...
                        F22'*reshape(F, [numBasisFun_Tar,numBasisFun_Src])*F11 - ...
                        F22*reshape(F, [numBasisFun_Tar,numBasisFun_Src])*F11' + ...
                        F22'*F22*reshape(F, [numBasisFun_Tar,numBasisFun_Src]),[],1),...
        F11_all, F22_all,'un',0)),2);
        
        if energy_settings.useLegacyNormalization
            alpha_nor = alpha/E_orientationTerm(C_eye);
            if ~isinf(alpha_nor)
                alpha = alpha_nor;
            end
        end
        
        E_orientationTerm = @(F) alpha*E_orientationTerm(F);
        dE_orientationTerm = @(F) alpha*dE_orientationTerm(F);
        
        Energies{end+1} = E_orientationTerm;
        dEnergies{end+1} = dE_orientationTerm;
    end

    if isempty(Energies)
        throw(MException('Energy:Null', 'No Energy was computed.'));
    end
    
end


function [ nor ] = compute_Res_nor(Ev1, Ev2, numEigsSrc, numEigsTar)
    % complex part
        x1 = Ev2./(1 + Ev2.^2);
        x2 = Ev1'./(1+ Ev1'.^2);
        Dlb1 = (repmat(x1, [1,numEigsSrc]) - repmat(x2, [numEigsTar,1])).^2;
        
        % real part
        x1 = 1./(1 + Ev2.^2);
        x2 = 1./(1+ Ev1'.^2);
        Dlb2 = (repmat(x1, [1,numEigsSrc]) - repmat(x2, [numEigsTar,1])).^2;
        % equal weight
        Dlb = Dlb1 + Dlb2;
        nor = norm(Dlb,'fro')
end

function [C,v,eval] = mat_projection(W)
n1 = size(W,1);
n2 = size(W,2);
[s,v,d] = svd(full(W));
C = s*eye(n1,n2)*d';
v = diag(v);
eval = [range(v),mean(v),var(v)];
end
