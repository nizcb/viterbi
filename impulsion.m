
d0 = 1;
d1 = 100;
N = 100; 


EbN0dB_min  = -2; 
EbN0dB_max  = 10; 
EbN0dB_step = 1;  
EbN0dB  = EbN0dB_min:EbN0dB_step:EbN0dB_max; 
trellis = poly2trellis(3,[7 5],7); 
Nombre_bits = 1024; 
numStates = trellis.numStates;
m = log2(numStates);
N= 2*(Nombre_bits+m);


longeur_EBN0 = length(EbN0dB);
TEP = zeros(1,longeur_EBN0);
for i=1:longeur_EBN0
   TEP(i) = methode_impulsion(d0, d1, N, Nombre_bits, trellis, EbN0dB(i));
end


figure;
semilogy(EbN0dB, TEP, 'LineWidth', 1.5);
xlabel('Eb/N0 (dB)');
ylabel('TEP');
title('TEP avec la méthode impulsion');
grid on; 


function TEP = methode_impulsion(d0, d1, N, Nombre_bits, trellis, EbN0dB)
    Eb_N0 = 10^(EbN0dB/10); 
    R = Nombre_bits/N; 
    
    
    y = ones(1, N); 
    v = zeros(1, Nombre_bits); 
    x_u = zeros(1, Nombre_bits); 
   
    for l = 1:Nombre_bits
        A = d0 - 0.5;
        x = x_u;
        while ((A <= d1) && isequal(x, x_u))
            A = A + 1;
            y(1 + (l-1)*2) = 1 - A; 
            x = viterbi_decode(y, trellis); 
        end
        v(l) = floor(A); 
    end

        D = unique(v);
    TEP = 0; 
    
    
    for d = D
        Ad = sum(v == d); 
        if d * R * Eb_N0 > 0
            TEP = TEP + (1/2) * Ad * erfc(sqrt(d * R * Eb_N0)); 
        else
            
            warning('L''argument de erfc est négatif ou zéro pour d = %d.', d);
        end
    end
end