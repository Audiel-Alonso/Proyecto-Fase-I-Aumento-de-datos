function [labels, centroids_rgb] = kmeans_lab(imagen_rgb, k, usar_solo_color)
    % imagen_rgb: Imagen de entrada en formato RGB
    % k: Número de grupos
    % usar_solo_color: (Booleano) Si es true, ignora el brillo (L)
    
    if nargin < 3, usar_solo_color = false; end

    % 1. Convertir al espacio de color L*a*b*
    % Este espacio separa la luminosidad (L) del color (a, b)
    img_lab = rgb2lab(imagen_rgb);
    [filas, cols, ~] = size(img_lab);
    
    % Reorganizar los datos
    datos = reshape(img_lab, [], 3);
    
    % Si el usuario quiere ignorar sombras/brillos, usamos solo canales 2 y 3
    if usar_solo_color
        puntos = datos(:, 2:3); 
    else
        puntos = datos;
    end
    
    num_puntos = size(puntos, 1);
    
    % 2. Inicialización aleatoria de centroides
    indices = randperm(num_puntos, k);
    centroids = puntos(indices, :);
    
    labels = zeros(num_puntos, 1);
    convergencia = false;
    
    % 3. Ciclo de K-means
    while ~convergencia
        old_labels = labels;
        
        % Distancia euclidiana (en L*a*b* esta distancia es perceptual)
        distancias = zeros(num_puntos, k);
        for i = 1:k
            diff = puntos - centroids(i, :);
            distancias(:, i) = sum(diff.^2, 2);
        end
        
        [~, labels] = min(distancias, [], 2);
        
        % Actualización de promedios
        for i = 1:k
            puntos_grupo = puntos(labels == i, :);
            if ~isempty(puntos_grupo)
                centroids(i, :) = mean(puntos_grupo, 1);
            end
        end
        
        if isequal(labels, old_labels)
            convergencia = true;
        end
    end
    
    % 4. Calcular los colores promedio finales en RGB para visualización
    % Sacamos el promedio de los pixeles originales (RGB) basado en las nuevas etiquetas
    pixeles_rgb = double(reshape(imagen_rgb, [], 3));
    centroids_rgb = zeros(k, 3);
    for i = 1:k
        centroids_rgb(i, :) = mean(pixeles_rgb(labels == i, :), 1);
    end
    
    % Devolver matriz de etiquetas con la forma original
    labels = reshape(labels, filas, cols);
end