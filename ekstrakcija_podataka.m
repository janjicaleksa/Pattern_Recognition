function [slovo_vx,slovo_vy,slovo_p,slovo_x_poz,slovo_y_poz] = ekstrakcija_podataka(slovo)

    slovo_vx = slovo(1,:);
    slovo_vy = slovo(2,:);
    slovo_p = slovo(3,:);
    slovo_x_poz = cumsum(slovo_vx);
    slovo_y_poz = cumsum(slovo_vy);
    
end

