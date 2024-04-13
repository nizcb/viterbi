function u = viterbi_decode(y, trellis)
    y(y>=0)=1;
    y(y<0)=-1;
    
    numOutputSymbols = trellis.numOutputSymbols;
    n_b = log2(numOutputSymbols); 
    Longeur = floor(length(y) / n_b);
    numStates = trellis.numStates;
    m = log2(numStates);
    etat_suivant = trellis.nextStates;
    outputs = trellis.outputs; 

    
    u = zeros(1, Longeur);


    metriques = ones(numStates, 1 + Longeur);
    inputs_matrice = ones(numStates, 1 + Longeur) ;
    matrice_etats_precedents = ones(numStates, 1 + Longeur);
    metriques(1, 1) = 0;

   
    for x = 1:Longeur
        for i = 1:numStates
                for z = 1:size(outputs, 2)
                    k = etat_suivant(i, z);
                    l = outputs(i, z);   
                    output = de2bi(l, 'left-msb', n_b)';
                    etat_courant = metriques(k + 1, x + 1);
                    metriques(k + 1, x + 1) = min(metriques(k+1, x+1), metriques(i, x) + y(1 + n_b * (x - 1):n_b * x) * output);
                    if etat_courant ~= metriques(k + 1, x + 1)
                        matrice_etats_precedents(k + 1, x + 1) = i;
                        inputs_matrice(k + 1, x + 1) = z - 1;  
                end
            end
        end
    end

    [~, alpha] = min(metriques(:, Longeur + 1));
   

    % Reconstitution 
    for i = 1:Longeur
        u(Longeur - i + 1) = inputs_matrice(alpha, 2 + Longeur - i);
        alpha = matrice_etats_precedents(alpha, 2 + Longeur - i);  
    end
    u=u(1:end-m);
end
