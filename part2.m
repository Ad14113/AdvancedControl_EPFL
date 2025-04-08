data = load('data_position.mat');
data = data.data;

Gfs = cell(1, 4);
Gs = cell(1, 4);

Ts = 0.02;

for i = 1:1:4

    experiment_name = data { i }. name ;
    y = data { i }. y ;
    u = data { i }. u ;
    
    Z = detrend ( iddata (y , u , Ts , "Period" , 8191) ) ;
    
    freqs = ( pi /4096: pi /4096: pi ) / Ts ;
    Gf = spa ( Z , 8191 , freqs );

    Gfs{i} = Gf;
    
    y_derivative = lsim (1 - tf ('z', Ts ) ^ -1 , y ) ;
    Z = detrend ( iddata ( y_derivative , u , Ts , "Period" , 8191) ) ;
    
    z = tf ('z', Ts ) ;
    G_derivative = oe (Z , [10 , 10 , 1]) ;
    G = G_derivative / (1 - z ^ -1) ;

    Gs{i} = G;

end

%for i=1:1:4
%    figure;
%    bode(Gfs{i});
%    experiment_name = sprintf('Experiment %d', i);
%    title(experiment_name);
%end


errors_Gfs = zeros(length(Gfs));
error_mean = zeros(1, length(Gfs));

Gfs_mean = (Gfs{1} + Gfs{2} + Gfs{3} + Gfs{4})/4;


for i=1:1:length(Gfs)
    for j=1:1:length(Gfs)
        result = norm(Gfs{j}/Gfs{i} - 1, Inf);
        errors_Gfs(i, j) = result;
    end
    result_mean = norm(Gfs{i}/Gfs_mean - 1, Inf);
    error_mean(i) = result_mean;
end

max_Gfs = zeros(1, length(Gfs));

for i=1:1:length(Gfs)
    [max_i, idx_i] = max(errors_Gfs(i, :));
    max_Gfs(i) = max_i;
end

[max_mean, idx_mean] = max(error_mean);

[min_max_Gfs, idx_min_max_Gfs] = min(max_Gfs);

if min_max_Gfs < max_mean
    idx_nom = idx_min_max_Gfs;
    Gnom = Gfs{idx_nom};
else
    Gnom = Gfs_mean;

end

order = 11;  % ordre souhaité
Gnom_ss = fitfrd(Gnom, order);   % produit un objet ss

figure;
h = bodeplot(Gnom_ss, Gnom);
setoptions(h, 'PhaseVisible', 'off');   

Gmm = stack(1, Gfs{1}, Gfs{2}, Gfs{3}, Gfs{4});
[Gu, info] = ucover(Gmm, Gnom_ss, 9, 'InputMult');


W2 = info.W1;
W2opt = info.W1opt;

% info.W1 (forced to respect the degree N)
% info.W1opt -> gives the best W2 (in terms of precision) but can be too complex 

figure;
h = bodeplot(W2, Gfs{1}/Gnom_ss - 1, Gfs{2}/Gnom_ss - 1, Gfs{3}/Gnom_ss - 1, Gfs{4}/Gnom_ss - 1);
setoptions(h, 'PhaseVisible', 'off');   


W1s = makeweight(10, 1, 0.51);
W1s_inverse = W1s^-1;


Ts = 0.02;                      % période d’échantillonnage
W1z = c2d(W1s, Ts, 'tustin');   % méthode bilinéaire (Tustin)

W3s = makeweight(1000, 0.1, 0.1);  % Gain élevé en BF, donc forte pénalité, le
% but est que U(jw) ait une faible norme en basse fréquence afin de
% compenser la grande norme de R(jw) en basse fréquence et ainsi avoir une
% faible commande 
Ts = 0.02;                      % période d’échantillonnage
W3z = c2d(W3s, Ts, 'tustin');   % méthode bilinéaire (Tustin)

% Maintenant on veut desing un H inf controller de sorte à minimiser la
% norme infinie à la fois de W1 S afin d'avoir une performance robuste et
% de W2 T afin d'avoir une stabilité robuste 
% on ajoute également un W3 pour le U et on le met à 1/5

Gnom_ss = ss(Gnom_ss);
W1_ss = ss(W1z);
W2_ss = ss(W2);
W3_ss = ss(W3z);

K=mixsyn(Gnom_ss,W1_ss,W3_ss,W2_ss);  % On ne fait rien pour le fonction de transfert U 
% On calcule tout  partir du Gnominal 

S = feedback(1, Gnom_ss*K);
T = feedback(Gnom_ss*K,1);
U = feedback(K,Gnom_ss);

figure;
bodemag(S, 'b', inv(W1_ss), 'k--');
legend('S', 'W1^{-1}');

figure;
bodemag(U, 'r', inv(W3_ss), 'k--');
legend('U', 'W3^{-1}');

figure;
bodemag(T, 'g', inv(W2_ss), 'k--');
legend('T', 'W2^{-1}');

t = 0:Ts:5;

figure;
step(T, t);
title('Réponse indicielle de T (boucle fermée : sortie)');

figure;
step(U, t);
title('Réponse indicielle de U (commande)');



perf_norm = norm(W1_ss * S, Inf);
stab_norm = norm(W2_ss * T, Inf);

disp('------------------------------------------');
disp(['Norme infinie de W1 * S (performance) : ', num2str(perf_norm)]);
disp(['Norme infinie de W2 * T (stabilité robuste) : ', num2str(stab_norm)]);
disp('------------------------------------------');