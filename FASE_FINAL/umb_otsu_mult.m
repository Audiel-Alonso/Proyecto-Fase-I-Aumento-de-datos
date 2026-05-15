function U = umb_otsu_mult(I, L)

    [M, N] = size(I);
    U = zeros(M, N);

    % Intensidades
    intensidad = 0:L-1;

    % Histograma
    [pi, ~] = hist(double(I(:)), intensidad);
    pi = pi / (M*N);

    % Sumas acumuladas
    P = cumsum(pi);                  % Probabilidad acumulada
    m = cumsum(intensidad .* pi);    % Medias acumulada

    % Media global
    m_G = m(end);

    % Inicialización
    maxSigma = -inf;
    K1_opt = 0;
    K2_opt = 0;

    % Búsqueda
    for K1 = 1:L-3

        % Valores para K2 válidos
        K2_range = K1+1:L-2;

        % Probabilidades
        P_1 = P(K1);
        P_2 = P(K2_range) - P(K1);
        P_3 = 1 - P(K2_range);

        % Medias
        m1 = m(K1) / P_1;
        m2 = (m(K2_range) - m(K1)) ./ P_2;
        m3 = (m_G - m(K2_range)) ./ P_3;

        % Varianza entre clases
        sigma = P_1*(m1 - m_G)^2 + P_2 .* (m2 - m_G).^2 + P_3 .* (m3 - m_G).^2;

        % Máximo K1 actual
        [sigma_max_local, idx] = max(sigma);

        if sigma_max_local > maxSigma
            maxSigma = sigma_max_local;
            K1_opt = K1;
            K2_opt = K2_range(idx);
        end
    end

    % Segmentación
    U(I <= K1_opt) = 2;
    U(I > K1_opt & I <= K2_opt) = 1;
    U(I > K2_opt) = 0;

end