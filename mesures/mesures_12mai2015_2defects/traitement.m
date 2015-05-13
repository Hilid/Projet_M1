clear all
close all
clc


defect_low = load('defect_low.txt');
defect_high = load('defect_high.txt');

no_defect_low = load('no_defect_low.txt');
no_defect_high = load('no_defect_high.txt');


defect_low_norm = defect_low(:,2) ./ max(defect_low(:,2));
defect_high_norm = defect_high(:,2) ./ max(defect_high(:,2));
no_defect_low_norm = no_defect_low(:,2) ./ max(no_defect_low(:,2))
no_defect_high_norm = no_defect_high(:,2) ./ max(no_defect_high(:,2));


figure(1)
plot(defect_high(:,1), defect_high_norm, 'o-r', 'Linewidth',3)
hold on
plot(defect_low(:,1), defect_low_norm,'Linewidth',3);
legend('defect with high lvl','defect with low lvl');

figure(2)
plot(no_defect_high(:,1), no_defect_high_norm, 'o-r','Linewidth',3)
hold on
plot(no_defect_low(:,1), no_defect_low_norm,'Linewidth',3);
legend('no_defect with high lvl','no_defect with low lvl');

figure(3)

plot(defect_low(:,1), defect_low_norm,'o-r','Linewidth',3);
hold on
plot(no_defect_low(:,1), no_defect_low_norm,'Linewidth',3);
legend('defect','no defect');
