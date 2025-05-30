function [rhoTest,rhoGen,POVMsATest,POVMsAGen,dimsA,dimsAPrime,probBlockCondTest,probBlockCondGen,eigenValueMat,deltaMat] = ...
    create46PhaseImperfectStatesandPOVMsAliceWithDelta(probTest,phaseRandomQuality,decoys,probDecoysCondTest,probDecoysCondGen,...
    probsSignalConTest,probsSignalConGen,photonCutoff,Na,tol)
    arguments
        probTest (1,1) double {mustBeInRange(probTest, 0, 1)}
        phaseRandomQuality (1,1) double {mustBeInRange(phaseRandomQuality, 0, 1)}
        decoys (:,1) double {mustBePositive}
        probDecoysCondTest (:,1) double {mustBeProbDist}
        probDecoysCondGen (:,1) double {mustBeProbDist}
        probsSignalConTest (:,1) double {mustBeProbDist}
        probsSignalConGen (:,1) double {mustBeProbDist}
        photonCutoff (1,1) double {mustBePositive,mustBeInteger}
        Na (1,1) double {mustBePositive,mustBeInteger} = 17;
        tol (1,1) double {mustBeNonnegative} = 1e-15;
    end

    %convert all inputs to vpa
    probTest = vpa(probTest);
    phaseRandomQuality = vpa(phaseRandomQuality);
    decoys = vpa(decoys);
    probDecoysCondTest = vpa(probDecoysCondTest);
    probDecoysCondGen = vpa(probDecoysCondGen);
    probsSignalConTest = vpa(probsSignalConTest);
    probsSignalConGen = vpa(probsSignalConGen);

    %We need to calculate all quantities until photonCutoff+1 blocks
    photonCutoffPlus = photonCutoff+1;

    %calculate states and eigen decompositions
    [rhoAprimeCell,eigenVectorCell,eigenValueCell] = ...
        createSignalStatesImperfectPhase(decoys,phaseRandomQuality,Na,tol);
    
    %perform Schmidt-decomposition and calculate Alice's POVMs (not conditioned
    %on gen/test yet)
    [POVMsA,signalStatesAAPrime,probBlock,eigenValueMat,Vcell] = ...
        POVMSandSchmidtDecompStates(probTest,probsSignalConTest,probsSignalConGen,...
        decoys,probDecoysCondTest,probDecoysCondGen,eigenVectorCell,eigenValueCell,photonCutoffPlus,tol);
    
    %calculate final test and gen states and POVMs
    [statesGen,statesTest,POVMsAGen,POVMsATest,dimsA,dimsAPrime,probBlockCondGen,probBlockCondTest] = ...
        createGenandTestStatesandPOVMS(POVMsA,signalStatesAAPrime,probTest,...
        decoys,probDecoysCondTest,probDecoysCondGen,probBlock,photonCutoffPlus,tol);

    %return states gen and test as density matrices
    rhoGen = cellfun(@(x) removeNumImprecHermMat(x*x',tol), statesGen, 'UniformOutput', false);
    rhoTest = cellfun(@(x) removeNumImprecHermMat(x*x',tol), statesTest, 'UniformOutput', false);
    
    %only from probBlockCondTest we need last entry and we store it in a
    %new variable (coresponds to p(photonCutoff+1 block|test)
    % probBlockCondTestNPlus = probBlockCondTest(photonCutoff+2);
    deltaMat = spectralGap(eigenValueMat);

    %We need to return all quantities until photonCutoff blocks, i.e. we
    %remove last entry in all variables
    rhoTest = rhoTest(1:photonCutoff+1);
    rhoGen = rhoGen(1:photonCutoff+1);
    POVMsATest = POVMsATest(:,1:photonCutoff+1);
    POVMsAGen = POVMsAGen(:,1:photonCutoff+1);
    dimsA = dimsA(1:photonCutoff+1);
    dimsAPrime = dimsAPrime(1:photonCutoff+1);
    probBlockCondTest = probBlockCondTest(1:photonCutoff+1);
    probBlockCondGen = probBlockCondGen(1:photonCutoff+1);
    eigenValueMat = eigenValueMat(:,1:photonCutoff+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions creating states, performing Schmidt-Decomposition etc

%%
%Calculate states (pure) conditioned on test and gen and POVM elements 
% conditioned on test and gen. Furthermore, finds the dimensions of A,
% APrime and the probabilities of each block conditioned on test and gen.
function [statesGen,statesTest,POVMsAGen,POVMsATest,dimsA,dimsAPrime,probBlockCondGen,probBlockCondTest] = ...
    createGenandTestStatesandPOVMS(POVMsA,signalStatesAAPrime,probTest,decoys,probDecoysCondTest,probDecoysCondGen,probBlock,photonCutoff,tol)
arguments
        POVMsA (:,:) cell {mustBeNonempty,mustBeCellOf(POVMsA,"double")}
        signalStatesAAPrime (1,:) cell {mustBeNonempty,mustBeCellOf(signalStatesAAPrime,"double")}
        probTest (1,1) double {mustBeInRange(probTest, 0, 1)}
        decoys (:,1) double {mustBePositive}
        probDecoysCondTest (:,1) double {mustBeProbDist}
        probDecoysCondGen (:,1) double {mustBeProbDist}
        probBlock (1,:) double {mustBeNonnegative,mustBeNonempty} 
        photonCutoff (1,1) double {mustBePositive,mustBeInteger} = 2;
        tol (1,1) double {mustBeNonnegative} = 1e-15;
    end

    %dimensions of A and APrime per block
    dimsA = arrayfun(@(x) size(POVMsA{1,x},1), 1:photonCutoff+1);
    dimsAPrime = arrayfun(@(x) size(signalStatesAAPrime{1,x},1), 1:photonCutoff+1)./dimsA;
    
    %number of decoy intensities
    numDecoy = numel(decoys);
    
    %Condition states and POVM elements on test and gen rounds
    indexGen = nonzeros(kron(logical(probDecoysCondGen),[ones(2,1);zeros(2,1)]).*(1:4*numDecoy)');
    indexTest = nonzeros(kron(logical(probDecoysCondTest),[zeros(2,1);ones(2,1)]).*(1:4*numDecoy)');
    
    %collect all POVM elements in gen and test round
    POVMsAGenUnCond = POVMsA(indexGen(:),:); %gen
    POVMsATestUnCond = POVMsA(indexTest(:),:); %test
    
    %% Generation rounds

    %projectors on gen round per block
    totalGenProj = cell(1,photonCutoff+1);
    totalGenProj(1,:) = cellfun(@(x)zeros(size(x)),POVMsA(1,:),'UniformOutput',false);

    for index = 1:size(POVMsAGenUnCond,1)
        for numPhoton = 0:photonCutoff
            totalGenProj(1,numPhoton+1) = {totalGenProj{1,numPhoton+1} + POVMsAGenUnCond{index,numPhoton+1}};
        end
    end

    %Pr(gen|block k)
    probGenConBlock = arrayfun(@(x) trace(kron(totalGenProj{x},eye(dimsAPrime(x)))*...
        signalStatesAAPrime{x}*signalStatesAAPrime{x}'), 1:1:photonCutoff+1 );

    %Diagonalize generation projectors
    %unitary TF and resulting diagonal matrices by diagonalizing gen projectors
    unTF = cell(1,photonCutoff+1);
    Dgen = cell(1,photonCutoff+1);

    for indexBlock = 1:photonCutoff+1
        [UMat,DMat] = eig(totalGenProj{indexBlock});

        %remove small numerical imprecision
        [DMat, DMatin] = removeNumImprecUnitary(DMat,tol);
        [UMat, UMatIn] = removeNumImprecUnitary(UMat,tol);
        
        %Store matrices
        unTF{indexBlock} = UMat;
        Dgen{indexBlock} = DMat;
    end

    %calculate sqrt of Dgen
    sqrtDgen = cellfun(@(x) diag(sqrt(diag(x))), Dgen,'UniformOutput',false);
    invSqrtDgen = cellfun(@(x) diag(safeInverse(diag(x))), sqrtDgen,'UniformOutput',false);
    
    %Alice POVMs in Gen rounds
    POVMsAGen = POVMsAGenUnCond;
    for indexBlock = 1:photonCutoff+1
        for indexA = 1:numel(indexGen)
            Mat = invSqrtDgen{indexBlock}'*unTF{indexBlock}'*POVMsAGenUnCond{indexA,indexBlock}*unTF{indexBlock}*invSqrtDgen{indexBlock};
    
            %remove small numerical imprecision
            [Mat, MatIn] = removeNumImprecHermMat(Mat,tol);

            %Store POVMs
            POVMsAGen{indexA,indexBlock} = Mat;
        end
    end

    % calculate states and blocks conditioned on gen rounds
    statesGen = arrayfun(@(x) kron(sqrtDgen{x}*unTF{x}',eye(dimsAPrime(x)))*signalStatesAAPrime{x},1:1:photonCutoff+1,'UniformOutput',false);
    statesGen = cellfun(@(x) 1/sqrt(trace(x*x'))*x,statesGen,'UniformOutput',false);

    probBlockCondGen = (probGenConBlock.*probBlock)./(1-probTest);
    
    %% Test rounds

    %projectors on gen round per block
    totalTestProj = cell(1,photonCutoff+1);
    totalTestProj(1,:) = cellfun(@(x)zeros(size(x)),POVMsA(1,:),'UniformOutput',false);
    
    for index = 1:size(POVMsATestUnCond,1)
        for numPhoton = 0:photonCutoff
            totalTestProj(1,numPhoton+1) = {totalTestProj{1,numPhoton+1} + POVMsATestUnCond{index,numPhoton+1}};
        end
    end
    
    %Pr(test|block k)
    probTestConBlock = arrayfun(@(x) trace(kron(totalTestProj{x},eye(dimsAPrime(x)))*...
        signalStatesAAPrime{x}*signalStatesAAPrime{x}'), 1:1:photonCutoff+1 );


    %Diagonalize testing projectors
    %unitary TF from gen rounds since [proj_genn,proj_test] = 0  
    % resulting diagonal matrices of testing projectors
    Dtest = cell(1,photonCutoff+1);

    for indexBlock = 1:photonCutoff+1
        DMat = unTF{indexBlock}'*totalTestProj{indexBlock}*unTF{indexBlock};

        %remove small numerical imprecision
        [DMat, DMatin] = removeNumImprecUnitary(DMat,tol);
        
        %Store matrices
        Dtest{indexBlock} = DMat;
    end

    %calculate sqrt of Dtest
    sqrtDtest = cellfun(@(x) diag(sqrt(diag(x))), Dtest,'UniformOutput',false);
    invSqrtDtest = cellfun(@(x) diag(safeInverse(diag(x))), sqrtDtest,'UniformOutput',false);

    %Alice POVMs in Test rounds
    POVMsATest = POVMsATestUnCond;
    for indexBlock = 1:photonCutoff+1
        for indexA = 1:numel(indexTest)
            Mat = invSqrtDtest{indexBlock}'*unTF{indexBlock}'*POVMsATestUnCond{indexA,indexBlock}*unTF{indexBlock}*invSqrtDtest{indexBlock};
            
            %remove small numerical imprecision
            [Mat, MatIn] = removeNumImprecHermMat(Mat,tol);
            
            %Store POVMs
            POVMsATest{indexA,indexBlock} = Mat;
        end
    end

    % calculate states and blocks conditioned on gen rounds
    statesTest = arrayfun(@(x) kron(sqrtDtest{x}*unTF{x}',eye(dimsAPrime(x)))*signalStatesAAPrime{x},1:1:photonCutoff+1,'UniformOutput',false);
    statesTest = cellfun(@(x) 1/sqrt(trace(x*x'))*x,statesTest,'UniformOutput',false);

    probBlockCondTest = (probTestConBlock.*probBlock)./probTest;
end


%% 
%Perform Schmidt-decomposition on each block and calculate POVM elements
function [POVMsA,signalStatesAAPrime,probBlock,lambdaMatk,Vcell] = POVMSandSchmidtDecompStates(probTest,probsSignalConTest,probsSignalConGen,decoys,probDecoysCondTest,probDecoysCondGen,eigenVectorCell,eigenValueCell,photonCutoff,tol)
    arguments
        probTest (1,1) double {mustBeInRange(probTest, 0, 1)}
        probsSignalConTest (:,1) double {mustBeProbDist}
        probsSignalConGen (:,1) double {mustBeProbDist}
        decoys (:,1) double {mustBePositive}
        probDecoysCondTest (:,1) double {mustBeProbDist}
        probDecoysCondGen (:,1) double {mustBeProbDist}
        eigenVectorCell cell {mustBeNonempty,mustBeCellOf(eigenVectorCell,"double")}
        eigenValueCell cell {mustBeNonempty,mustBeCellOf(eigenValueCell,"double")}
        photonCutoff (1,1) double {mustBePositive,mustBeInteger} = 2;
        tol (1,1) double {mustBeNonnegative} = 1e-15;  
    end

    %number of decoy intensities
    numDecoy = numel(decoys);

    %Create empty cell storing the transformations on system APrime
    Vcell = cell(1,photonCutoff+1);
    
    %Create empty cell for Alice's POVm elements (colums index block k)
    POVMsA = cell(4*numDecoy,photonCutoff+1);

    %Create empty cell for signalstates in AAPrime (colums index block k)
    signalStatesAAPrime = cell(1,photonCutoff+1);

    %create vector of zeros for probabilities of each block
    probBlock = zeros(photonCutoff+1,1);

    %probabilities of sending signals in test and gen rounds padded with 0s
    probsSignalGenExtended = (1-probTest)*[probsSignalConGen;0;0];
    probsSignalTestExtended = probTest*[0;0;probsSignalConTest];
    
    %probability of signal choice and intensity
    probSignalandIntensity = zeros(4*numDecoy,1);

    for indexIntensity = 1:numDecoy
        for indexStateA = 1:4
            probSignalandIntensity((indexIntensity-1)*4+indexStateA) = ...
                probsSignalTestExtended(indexStateA)*probDecoysCondTest(indexIntensity)+...
                    probsSignalGenExtended(indexStateA)*probDecoysCondGen(indexIntensity);
        end
    end
    
    %eigenvalues for each signal and intensity conditioned on k
    lambdaMatk = zeros(4*numDecoy,photonCutoff+1);
    
    for numPhoton = 0:photonCutoff
        temp = cell2mat(cellfun(@(x) x(numPhoton+1), eigenValueCell,'UniformOutput',false));
        lambdaMatk(:,numPhoton+1) = reshape(temp,[],1);
    end
    
    %prob
    omega = lambdaMatk.'*probSignalandIntensity;

    for numPhoton = 0:photonCutoff
        %Calculate unnormalized signal state
        unNormalizedSignalState = 0;
        for indexIntensity = 1:numDecoy
            for indexStateA = 1:4            
                %select eigenvectors corresponding to signalA and intensity
                Wk = eigenVectorCell{indexStateA,indexIntensity};
                %select eigenvalues corresponding to signalA and intensity
                lambdak = eigenValueCell{indexStateA,indexIntensity};
                
                %usual signal and intensity choice probabilities                
                prob = probSignalandIntensity((indexIntensity-1)*4+indexStateA)*...
                    lambdaMatk((indexIntensity-1)*4+indexStateA,numPhoton+1);
        
                %Construct normalized signal state
                stateA = zket(4*numDecoy,(indexIntensity-1)*4+indexStateA);
                stateAPrime = Wk(:,numPhoton+1);
                unNormalizedSignalState = unNormalizedSignalState +...
                    sqrt(prob/omega(numPhoton+1))*kron(stateA, stateAPrime);
            end
        end
        
        %probability of block k
        probk = trace(unNormalizedSignalState*unNormalizedSignalState');
        probBlock(numPhoton+1) = omega(numPhoton+1);

        %normalized signal state
        normalizedSignalState = 1/sqrt(probk)*unNormalizedSignalState;
        
        %Perform Shmidt-decomposition
        [singVals,U,V] = SchmidtDecomposition(normalizedSignalState,[4*numDecoy,size(Wk(:,numPhoton+1),1)]);
        
        %If errors occur introduce cutoff for singular values!!!!

        %remove imprecision from U and V
        [U,UIn] = removeNumImprecUnitary(U,tol);
        [V,VIn] = removeNumImprecUnitary(V,tol);
        
        %store V for debugging
        Vcell{numPhoton+1} = V;

        %schmidt-rank
        srank = numel(singVals);
        
        %Schmidt reduced state in A and APrime
        signalStatesAAPrime{numPhoton+1} = kron(diag(singVals),eye(srank))*MaxEntangled(srank,0,0);
        
        %Create Alice's POVM elements
        for indexA = 1:4*numDecoy
            Mat = U'*diag(zket(4*numDecoy,indexA))*U;
            %remove small numerical imprecision
            [Mat, MatIn] = removeNumImprecHermMat(Mat,tol);
            POVMsA{indexA,numPhoton+1} = Mat;
        end
    end
end

%%
 %Creates the signal states with phase imperfection and performs the
    %eigendecomposition
function [rhoAprimeCell,eigenVectorCell,eigenValueCell] = ...
    createSignalStatesImperfectPhase(decoys,phaseRandomQuality,Na,tol)

    arguments
        decoys (:,1) double {mustBePositive}
        phaseRandomQuality (1,1) double {mustBeInRange(phaseRandomQuality, 0, 1)}
        Na (1,1) double {mustBePositive,mustBeInteger} = 17;
        tol (1,1) double {mustBeNonnegative} = 1e-15;  
    end

    %Creates the signal states with phase imperfection and performs the
    %eigendecomposition

    %convert decoys to vps
    decoys = vpa(decoys);
    
    %Number of decoy intensities
    numDecoy = numel(decoys);
    
    %Store signal states and eigenvectors and -values in cells for each
    %signal state choice and intensity
    rhoAprimeCell = cell(4,numDecoy);
    eigenVectorCell = cell(4,numDecoy);
    eigenValueCell = cell(4,numDecoy);
    
    for indexIntensityA = 1:numDecoy       
        %% HV Basis

        %complex amplitudes of cpoherent states
        complexAmplHV0 = [sqrt(decoys(indexIntensityA));0];
        complexAmplHV1 = [0;sqrt(decoys(indexIntensityA))];

        %0 bit
        cohMatHV0 = coherentStatePhotonNumBasis(complexAmplHV0,Na);
        fullPhaseMatHV0 = fullyPhaseRandomNumBasis(norm(complexAmplHV0)^2,0,0,Na);

        rhoAprimeCell{1,indexIntensityA} = ...
            phaseRandomQuality*fullPhaseMatHV0 + (1-phaseRandomQuality)*cohMatHV0;

        %calculate eigenvalues and -vectors
        [eigenVectorCell{1,indexIntensityA},eigenValueCell{1,indexIntensityA}] = ...
            sortedEigs(rhoAprimeCell{1,indexIntensityA},tol);
        
        %1 bit
        cohMatHV1 = coherentStatePhotonNumBasis(complexAmplHV1,Na);
        fullPhaseMatHV1 = fullyPhaseRandomNumBasis(norm(complexAmplHV1)^2,1,0,Na);

        rhoAprimeCell{2,indexIntensityA} = ...
            phaseRandomQuality*fullPhaseMatHV1 + (1-phaseRandomQuality)*cohMatHV1;

        %calculate eigenvalues and -vectors
        [eigenVectorCell{2,indexIntensityA},eigenValueCell{2,indexIntensityA}] = ...
            sortedEigs(rhoAprimeCell{2,indexIntensityA},tol);

        %% DA Basis
        %complex amplitudes of cpoherent states
        uTransDA = 1/sqrt(2)*[1,1;1,-1];

        complexAmplDA0 = uTransDA*complexAmplHV0;
        complexAmplDA1 = uTransDA*complexAmplHV1;

        %0 bit
        cohMatDA0 = coherentStatePhotonNumBasis(complexAmplDA0,Na);
        fullPhaseMatDA0 = fullyPhaseRandomNumBasis(norm(complexAmplDA0)^2,0,1,Na);

        rhoAprimeCell{3,indexIntensityA} = ...
            phaseRandomQuality*fullPhaseMatDA0 + (1-phaseRandomQuality)*cohMatDA0;

        %calculate eigenvalues and -vectors
        [eigenVectorCell{3,indexIntensityA},eigenValueCell{3,indexIntensityA}] = ...
            sortedEigs(rhoAprimeCell{3,indexIntensityA},tol);

        %1 bit
        cohMatDA1 = coherentStatePhotonNumBasis(complexAmplDA1,Na);
        fullPhaseMatDA1 = fullyPhaseRandomNumBasis(norm(complexAmplDA1)^2,1,1,Na);

        rhoAprimeCell{4,indexIntensityA} = ...
            phaseRandomQuality*fullPhaseMatDA1 + (1-phaseRandomQuality)*cohMatDA1;

        %calculate eigenvalues and -vectors
        [eigenVectorCell{4,indexIntensityA},eigenValueCell{4,indexIntensityA}] = ...
            sortedEigs(rhoAprimeCell{4,indexIntensityA},tol);
    end
end
%%
%Creates density matrix of coherent state in HV basis with <= Na photon cutoff
function cohMat = coherentStatePhotonNumBasis(alphaHV,Na)
    arguments
        alphaHV (1,2) double
        Na (1,1) double {mustBePositive,mustBeInteger} = 17;
    end
    %Creates density matrix of coherent state in HV basis with <= Na photon cutoff
    
    %Recast components of vector alphaHV
    alphaH = alphaHV(1);    
    alphaV = alphaHV(2);

    cohVecFull = 0;
    for numPh = 0:Na
        for numH = 0:numPh
            numV = numPh - numH;
            %coeff of |numH,numV> for coherent state |alphaH,alphaV>
            coeff = double(coherentStatePhotonDist(alphaH,numH)*coherentStatePhotonDist(alphaV,numV));
            
            %create vector in HV basis
            vector = converttoZ(Na,numH,numV,0);
    
            %add to coherent state vector
            cohVecFull = cohVecFull + coeff*vector;
        end
    end
    %Project on <= Na photon space
    projector = projector_Nph_subspace(Na);
    cohMat = projector*(cohVecFull*cohVecFull')*transpose(projector);
end

%helper to create photon-number distribution of a coherent state
function prob = coherentStatePhotonDist(alpha,n)
    arguments
        alpha (1,1) double
        n double {mustBeNonnegative,mustBeInteger}
    end
    %absolute value of alpha
    absAlphaSq = vpa(abs(alpha)^2);

    prob = exp(-absAlphaSq/2)*alpha^n/sqrt(factorial(n));
end

%%
%Creates density matrix of a fully phase-randomized state in HV basis with <= Na photon cutoff
function stateMat = fullyPhaseRandomNumBasis(intensity,value,basis,Na)
    arguments
        intensity (1,1) double {mustBeNonnegative}
        value (1,1) double {mustBeMember(value,[0,1])}
        basis (1,1) double {mustBeMember(basis,[0,1,2])}
        Na (1,1) double {mustBePositive,mustBeInteger} = 17;
    end
    %Creates density matrix of a fully phase-randomized state in HV basis with <= Na photon cutoff
    stateMatFull = 0;
    if value == 0
        for num = 0:Na
            %coeff of |numH,numV> for fully phase-randomized state
            coeff = coherentStatePhotonDist(sqrt(vpa(intensity)),num);

            %create vector in HV basis
            vector = converttoZ(Na,num,0,basis);
            stateMatFull = stateMatFull + double(abs(coeff)^2)* (vector*vector');
        end
    elseif value == 1
        for num = 0:Na
            %coeff of |numH,numV> for fully phase-randomized state
            coeff = coherentStatePhotonDist(sqrt(vpa(intensity)),num);
            
            %create vector in HV basis
            vector = converttoZ(Na,0,num,basis);
            stateMatFull = stateMatFull + double(abs(coeff)^2) * (vector*vector');
        end
    end 

    %Project on <= Na photon space
    projector = projector_Nph_subspace(Na);
    stateMat = projector*stateMatFull*transpose(projector);
end

%% 
%Convert state |n,m>_basis to z-basis with cutoff Nmax

function vector = converttoZ(Nmax,m,n,basis)
    arguments
        Nmax (1,1) double {mustBeNonnegative,mustBeInteger}
        m (1,1) double {mustBeNonnegative,mustBeInteger}
        n (1,1) double {mustBeNonnegative,mustBeInteger}
        basis (1,1) double {mustBeMember(basis,[0,1,2])}        
    end
    %takes a ket |m, n> in the X,Y,Z basis and converts to Z basis.
    vector = zeros((Nmax+1)^2,1);
    
    %Z basis
    if basis == 0
        Hspace = zeros(Nmax+1,1);
        Vspace = zeros(Nmax+1,1);

        Hspace(m+1) = 1;
        Vspace(n+1) = 1;

        vector = kron(Vspace, Hspace);
    
    %X basis
    elseif basis == 1        
        for l=0:m
            for k=0:n             
                coeffX = nchoosek(m,l) * nchoosek(n,k) * (-1)^(n-k) * sqrt(factorial(l+k) * factorial(m+n-l-k)/(factorial(n) *factorial(m)));

                Hspace = zeros(Nmax+1,1);
                if l+k <= Nmax                                                         
                    Hspace(l+k+1) = 1;                                       
                end

                Vspace = zeros(Nmax+1,1);
                if n+m-l-k <= Nmax
                    Vspace(m+n-l-k+1) = 1;
                end

                vector = vector + coeffX/sqrt(2)^(m+n) * kron(Hspace, Vspace);
            end
        end

    %Y basis
    elseif basis == 2
        for l=0:m
            for k=0:n
                coeffY = nchoosek(m,l) * nchoosek(n,k) * (-1)^(n-k) * (1i)^(m+n-l-k) * sqrt(factorial(l+k) * factorial(m+n-l-k)/(factorial(n) *factorial(m)));

                Hspace = zeros(Nmax+1,1);
                if l+k <= Nmax                                                         
                    Hspace(l+k+1) = 1;                                       
                end

                Vspace = zeros(Nmax+1,1);
                if n+m-l-k <= Nmax
                    Vspace(m+n-l-k+1) = 1;
                end

                vector = vector + coeffY/sqrt(2)^(m+n) * kron(Hspace, Vspace);
            end
        end
    end
end

%%
%Creates projector onto <= Na photon subspace

function projector = projector_Nph_subspace(Na)
    arguments
        Na (1,1) double {mustBePositive,mustBeInteger}
    end
    %Creates a projector onto to the <= Nb photon subspace
    %Nb: photon number cutoff

    n = 1; %number of spaces here n=1
       
    dim = (Na+1) * nchoosek(2*n+Na, 2*n-1)/(2*n); % sum (2*n+i-1) choose (2*n-1) from i = 0 to i = N
    
    projector = zeros(dim, (Na+1)^(2*n)); % Projector onto the <= N photon subspace with block diagonal
    projector(1,1) = 1;                  % structure corresponding to total photon number in n-bit block.
    
                                  
    for photonnumber = 1:1:Na % Considering all photon numbers till the photon number cut-off
        %% Forming vector of all the ways that $x_1+....+x_{2n} = photonnumber$
        temp = sym2cell(feval(symengine, 'combinat::compositions', photonnumber, strcat('Length = ',int2str(2*n)), 'MinPart = 0'));
        
        l = size(temp,2); % Count number of such combinations that are there.
            for i = 1:1:l % Loop over all combinations
                combtemp = sym2cell(temp{i});
                index = 1; % Column that will be non-zero in our projection matrix
                for j = 1:1:2*n % Loop over all the elements in the combination
                    index = index + combtemp{j}*(Na+1)^(j-1); % |....,i,j,k> has only the 1+(N+1)*k+(N+1)^2*j+(N+1)^3*i+... element non-zero
                                                            % Thus, this is the column of the projector that will be non-zero
                end
                projector((photonnumber) * nchoosek(2*n+photonnumber-1, 2*n-1)/(2*n)+i, index) = 1;
            end
    end
end

%% 
%Function sorts eigenvalues from largest to smallest
function [V,d] = sortedEigs(Mat,tol)
    arguments
        Mat (:,:) double 
        tol (1,1) double {mustBeNonnegative} = 1e-15;
    end
    %calculates sorted eigenvectors and -values by largest eigenvalue

    %Perform eigen decomposition (unsorted)
    [V,D] = eig(Mat);
    d = double(diag(D));
    
    %Sort eigenvalues and vectors in descending order
    if ~issorted(d,'descend')
        [d,inds] = sort(d,'descend');
        d(abs(d) <= tol ) = 0;
        V = V(:, inds);
        V(abs(V) <= tol ) = 0;
    end

    V = double(V);
end

%% 
%Function removes numerical imprecision from matrix and makes it hermitian,
%also warns if matrix is very non-hermitian

function [MatOut,MatIn] = removeNumImprecHermMat(MatIn,tol)
    arguments
        MatIn (:,:) double 
        tol (1,1) double {mustBeNonnegative} = 1e-15;
    end
    %Calculate maximum nonhermitian part
    maxNonHerm = max(abs(MatIn - MatIn'),[],"all");
    
    if maxNonHerm > 1e-6
        throwAsCaller(MException("removeNumImprec:MatrixNotHermitian",...
            "Matrix is not Hermitian."))
    end

    %Set output equal input
    MatOut = 1/2*(MatIn+MatIn');

    %set small values = 0
    MatOut(abs(MatOut) <= tol ) = 0;
    
    %set small imaginary values = 0
    imMat = imag(MatOut);
    imMat(abs(imMat) <= tol ) = 0;

    %set small real values = 0
    reMat = real(MatOut);
    reMat(abs(reMat) <= tol ) = 0;
    
    %put matrix back together
    MatOut = reMat + 1i*imMat;

    %guarantee MatOut is hermitian
    MatOut = 1/2*(MatOut + MatOut');
end

function [MatOut,MatIn] = removeNumImprecUnitary(MatIn,tol)
    arguments
        MatIn (:,:) double 
        tol (1,1) double {mustBeNonnegative} = 1e-15;
    end

    %Set output equal input
    MatOut = MatIn;

    %set small values = 0
    MatOut(abs(MatOut) <= tol ) = 0;
    
    %set small imaginary values = 0
    imMat = imag(MatOut);
    imMat(abs(imMat) <= tol ) = 0;

    %set small real values = 0
    reMat = real(MatOut);
    reMat(abs(reMat) <= tol ) = 0;
    
    %put matrix back together
    MatOut = reMat + 1i*imMat;
end

%calculates component-wise inverse on support of vector
function vecOut = safeInverse(vecIn)
    %assign all zeros to output
    vecOut = zeros(size(vecIn));
    %find non-zero elements in input
    nonZeroElmnts = nonzeros(vecIn);
    %find indices of non-zero elements in input
    indexNonZero = find(vecIn);
    %calc inverse of non-zero elements
    invNonZeroElmnts = 1./nonZeroElmnts;
    %assign inverse to non-zero elements 
    vecOut(indexNonZero) = invNonZeroElmnts;
end

%Poisson distribution
function prob = poisson(intensity, count)
    prob = exp(-intensity).*intensity.^count./factorial(count);
end

%calculate delta_m from eq.(207)
% probBlockCondTest needs to be ordered 0,1,... !!!!
function delta = spectralGap(probBlockCondAandMu)
    delta = zeros(size(probBlockCondAandMu,1),size(probBlockCondAandMu,2)-1);
    delta(:,1) = abs(probBlockCondAandMu(:,2) - probBlockCondAandMu(:,1));

    for index = 2:size(probBlockCondAandMu,2)-1
        delta(:,index) = min(abs(probBlockCondAandMu(:,index+1) - probBlockCondAandMu(:,index)), ...
            abs(probBlockCondAandMu(:,index) - probBlockCondAandMu(:,index-1)));
    end

    %make sure delta >=0 
    delta = max(delta,zeros(size(delta)));
end