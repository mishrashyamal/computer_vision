function seam = find_vertical_seam(energy_map)
    seam_mat = energy_map;

    for i = 2:size(seam_mat, 1)
        for j = 1:size(seam_mat, 2)
            if j == 1
                seam_mat(i, j) = seam_mat(i, j) + min(seam_mat(i-1, j), seam_mat(i-1, j+1));
            elseif j == size(seam_mat, 2)
                seam_mat(i, j) = seam_mat(i, j) + min(seam_mat(i-1, j-1), seam_mat(i-1, j));
            else
                seam_mat(i, j) = seam_mat(i, j) + min(seam_mat(i-1, j-1), min(seam_mat(i-1, j), seam_mat(i-1, j+1)));
            end
        end
    end

   
    [~, idx] = min(seam_mat(end, :));
    opt_seam = zeros(size(seam_mat, 1), 1);
    opt_seam(end) = idx;

    for i = size(seam_mat, 1)-1:-1:1
        if opt_seam(i+1) == 1
            [~, m] = min([seam_mat(i, 1), seam_mat(i, 2)]);
            opt_seam(i) = m;
        elseif opt_seam(i+1) == size(seam_mat, 2)
            [~, m] = min([seam_mat(i, end-1), seam_mat(i, end)]);
            opt_seam(i) = size(seam_mat, 2) - (m-1);
        else
            [~, m] = min([seam_mat(i, opt_seam(i+1)-1), seam_mat(i, opt_seam(i+1)), seam_mat(i, opt_seam(i+1)+1)]);
            opt_seam(i) = opt_seam(i+1) + (m-2);
        end
    end

    seam = opt_seam;
end
