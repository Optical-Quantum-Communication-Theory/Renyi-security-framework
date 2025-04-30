function POVMs = Create_POVM_6State_Flag(Nb,prob_list)
    % Creates cell array of Bob's POVM elements in the following order:
    % single clicks: H, V, +, -, R, L
    % multi clicks: mult
    % vacuum

    %Flag vectors follow convention:  flag H, flag V, flag +, flag -,
    %flag R, flag L, flag multiclick
    
    % Nb: photon number cut-off
    % prob_list: list of Bob's basis probabilities

    %create empty array for storing POVMs
    POVMs = cell(1,2*numel(prob_list)+2);
    
    %Projector onto <= Nb photon subspace
    projector = projector_Nph_subspace(Nb);

    %% Add single-clicks to POVM list in the order H, V, +, -, R, L
    for basis=1:numel(prob_list) 
        for value=1:2
            %Single-click POVM element in <= Nb subspace
            Fsingle = projector*single_click(Nb,prob_list(basis),value-1,basis-1)*transpose(projector);

            %Vector needed for flag space
            flagvecsingle = zket(2*numel(prob_list)+1, 2*(basis-1)+(value-1)+1);

            %Final POVM element including flag
            POVMs{2*(basis-1)+(value-1)+1} = 1/2*(blkdiag(Fsingle,diag(flagvecsingle)) + blkdiag(Fsingle,diag(flagvecsingle))');
        end
    end
    
    %% Add multiclick by claculating double and cross-clicks first
    %Double-clicks in the order HV, +-, RL
    Fdouble = cell(1,3);
    for basis = 1:numel(prob_list)
        %Double-click POVM element in <= Nb subspace
        FdoubleInit = projector*double_click(Nb,prob_list(basis),basis-1)*transpose(projector);
        
        %Final POVM element excluding flag
        Fdouble{basis} = 1/2*(FdoubleInit + FdoubleInit');
    end

    %Cross-clicks
    %Cross-click POVM element in <= Nb subspace
    FcrossInit = projector*cross_click(Nb,prob_list)*transpose(projector);
    Fcross = 1/2*(FcrossInit + FcrossInit');
    
    %Vector needed for flag space of MultiClick
    flagvecMult = zket(2*numel(prob_list)+1, 2*numel(prob_list)+1);

    %Final multi-click POVM element including flag
    Fmult = sum(cat(3,Fdouble{:}),3) + Fcross;
    Fmult = 1/2*(Fmult + Fmult');
    POVMs{2*numel(prob_list)+1} = blkdiag(Fmult,diag(flagvecMult));

    %% Add vacuum POVM element
    Fvac = zeros(size(POVMs{1}));
    Fvac(1,1) = 1;
    POVMs{2*numel(prob_list)+2} = Fvac;
end

function povm = single_click(Nb,basis_prob,value,basis)
    % Creates a single-click POVM element in particular basis with value
    % 0,1
    % Nb: photon number cut-off
    % basis_prob: basis probability
    % value: 0,1 e.g in Z basis H = 0, V = 1
    % basis: basis choice, 0 = Z, 1 = X, 2 = Y 

    povm = 0;
    if value ==0
        for num=1:Nb
            vector = converttoZ(Nb,num,0,basis);
            povm = povm + basis_prob^num * (vector*vector');
        end
    elseif value==1
        for num=1:Nb
            vector = converttoZ(Nb,0,num,basis);
            povm = povm + basis_prob^num * (vector*vector');
        end
    end   
end

function povm = double_click(Nb,basis_prob,basis)
    % Creates a double-click POVM element in particular basis
    % Nb: photon number cut-off
    % basis_prob: basis probability
    % basis: basis choice, 0 = Z, 1 = X, 2 = Y 

    povm = 0;
    for num = 2:Nb
        for k = 1:num-1
            vector = converttoZ(Nb,k,num-k,basis);
            povm = povm + basis_prob*(vector*vector');
        end
    end
end

function povm = cross_click(Nb,basis_prob_list)
    % Creates the cross-click POVM element
    % Nb: photon number cut-off
    % basis_prob_list: list basis probabilities

    povm = 0;
    for basis_elmnt=0:numel(basis_prob_list)-1
        for num= 1:Nb
            vector0 = converttoZ(Nb,num,0,basis_elmnt);
            vector1 = converttoZ(Nb,0,num,basis_elmnt);
            povm = povm + basis_prob_list(basis_elmnt+1)*(1-basis_prob_list(basis_elmnt+1)^(num-1)) * ( (vector0*vector0') + (vector1*vector1'));
        end
    end
end

function vector = converttoZ(Nmax,m,n,basis)
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

function projector = projector_Nph_subspace(Nb)
    %Creates a projector onto to the N <= Nb photon subspace
    %Nb: photon number cutoff

    n = 1; %number of spaces here n=1
       
    dim = (Nb+1) * nchoosek(2*n+Nb, 2*n-1)/(2*n); % sum (2*n+i-1) choose (2*n-1) from i = 0 to i = N
    
    projector = zeros(dim, (Nb+1)^(2*n)); % Projector onto the <= N photon subspace with block diagonal
    projector(1,1) = 1;                  % structure corresponding to total photon number in n-bit block.
    
                                  
    for photonnumber = 1:1:Nb % Considering all photon numbers till the photon number cut-off
        %% Forming vector of all the ways that $x_1+....+x_{2n} = photonnumber$
        temp = sym2cell(feval(symengine, 'combinat::compositions', photonnumber, strcat('Length = ',int2str(2*n)), 'MinPart = 0'));
        
        l = size(temp,2); % Count number of such combinations that are there.
            for i = 1:1:l % Loop over all combinations
                combtemp = sym2cell(temp{i});
                index = 1; % Column that will be non-zero in our projection matrix
                for j = 1:1:2*n % Loop over all the elements in the combination
                    index = index + combtemp{j}*(Nb+1)^(j-1); % |....,i,j,k> has only the 1+(N+1)*k+(N+1)^2*j+(N+1)^3*i+... element non-zero
                                                            % Thus, this is the column of the projector that will be non-zero
                end
                projector((photonnumber) * nchoosek(2*n+photonnumber-1, 2*n-1)/(2*n)+i, index) = 1;
            end
    end
end