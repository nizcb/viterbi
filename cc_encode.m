
function c = cc_encode(u, trellis)
    n_b = log2(trellis.numOutputSymbols);
    c = zeros(1, length(u) * n_b);
    etat_actuel = 0;
    etat_suivant = trellis.nextStates;
    sortie = trellis.outputs;
    numStates = trellis.numStates;
    m=log2(numStates);

    for i = 1:length(u)
        Symbol_entree = u(i) + 1;
        binaire = fliplr(de2bi(sortie(etat_actuel + 1, Symbol_entree), n_b));
        for j = 1:n_b
            c((i - 1) * n_b + j) = 1-2*binaire(j);
        end
        etat_actuel = etat_suivant(etat_actuel + 1, Symbol_entree);
    end
    for i = length(u) + 1 : length(u) + m
        if trellis.nextStates(etat_actuel+1,1) == floor(etat_actuel/2)
        binaire = fliplr(de2bi(sortie(etat_actuel + 1, 1), n_b));
        for j = 1:n_b
            c((i - 1) * n_b + j) = 1-2*binaire(j);
        end

        
         
        else
            binaire = fliplr(de2bi(sortie(etat_actuel + 1, 2), n_b));
        for j = 1:n_b
            c((i - 1) * n_b + j) = 1-2*binaire(j);
        end
        
        end
    end
end
