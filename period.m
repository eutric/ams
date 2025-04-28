function T = period(O)
% Funzione che dall'oggetto orbita in ingresso, restituisce il valore in
% secondi del periodo di questa e stampa il risultato in ore
T = 2*pi*sqrt(O.a^3/O.mu);
minuti = rem(T,3600);
fprintf('hh:mm - %d:%d',floor(T/3600),minuti/60)
end

