function [Ob1,Ob2,Ob3,Ob4,Ob5,Ob6,Ob7] = ekstrakcija_obelezja(slovo)

    [slovo_vx,slovo_vy,slovo_p,slovo_x_poz,slovo_y_poz] = ekstrakcija_podataka(slovo);

    Ob1 = mean(slovo_vx); % srednja brzina po x osi
    Ob2 = mean(slovo_vy); % srednja brzina po y osi
    Ob3 = max(slovo_x_poz) - min(slovo_x_poz); % sirina
    Ob4 = Ob1/(max(slovo_y_poz)-min(slovo_y_poz)); % odnos sirine i visine
    %broj odbiraka u donjoj polovini slike pola slike
    n1start = min(slovo_x_poz);
    n2start = min(slovo_y_poz);
    n1stop = max(slovo_x_poz);
    n2stop = max(slovo_y_poz)/2;
    Ob5 = histcounts2(slovo_x_poz,slovo_y_poz,[n1start n1stop],[n2start n2stop]);
    Ob6 = slovo_y_poz(1,1) - slovo_y_poz(1,length(slovo_y_poz)); %razlika pocetnog i krajnjeg polozaja po y osi
    Ob7 = max(slovo_p) - min(slovo_p); % razlika maksimalnog i minimalnog pritiska
    
end