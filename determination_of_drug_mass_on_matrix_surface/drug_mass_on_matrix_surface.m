clc
close all
clear all
pkg load image
REM=imread('REM_ohne_Legende.png');

% grayscale
REM_gray=rgb2gray(REM);
figure,
imshow(REM_gray);

% bw
threshold=0.78;
REM_bw=im2bw(REM_gray, threshold);
figure,
imshow(REM_bw);

% count
anzahl = bwconncomp(REM_bw);
anzahlPartikel = anzahl.NumObjects;
disp('Die Anzahl der einzelnen Medikamentenpartikel betraegt');
disp(anzahlPartikel)
disp('');
disp('');


disp('R E M  -  A U F N A H M E');

% area Partikel
area = bwarea(REM_bw);
disp('Die Flaeche, die alle auf dem Ausschnitt zu sehenden Partikel zusammen einnehmen, betraegt in Pixeln');
disp(area);
display('');

% area REM
area_REM_Pixel = 1072*730;
area_REM_mu2 = 221.58;
disp('Die REM-Aufnahme hat eine Auflösung von 1072x730 Pixeln und besteht damit aus ... Pixeln');
disp(area_REM_Pixel)
display('');
disp('Die Fläche, die alle auf dem Ausschnitt zu sehenden Partikel zusammen einnehmen, betraegt in mu^2');
area_mu2 = (area_REM_mu2/area_REM_Pixel)*area;
disp(area_mu2)
disp('');
disp('');


disp('P A R T I K E L G E O M E T R I E');

# durchschnittlicher Partikelquerschnitt
area_mu2_pro_partikel = area_mu2/anzahlPartikel;
disp('Der durchschnittliche Partikelquerschnitt beträgt ...[mu m²].');
disp(area_mu2_pro_partikel)
disp('');

% area Partikel Prozent
area_prozent = area/area_REM_Pixel*100;
disp('Damit nehmen die Medikamentenpartikel ... Prozent der Oberflaeche ein.');
disp(area_prozent)
disp('');

# durchschnittlicher Partikelradius
disp('Der durchschnittliche Partikelradius beträgt ...[mu m].');
r = sqrt(area_mu2_pro_partikel/pi);
disp(r)
disp('');

# durchschnittliche Partikeloberfläche
disp('Die durchschnittliche Partikeloberfläche beträgt ...[mu m^2].');
a = 4*pi*r^2;
disp(a)
disp('');
disp('');


disp('M E D I K A M E N T E N L A D U N G');

# Medikamentenvolumen auf REM-Aufnahme
disp('Das kummulierte Medikamentenvolumen auf der REM-Aufnahme beträgt ...[mu m^3].');
V_REM = (4/3)*pi*r^3*anzahlPartikel;
disp(V_REM)
disp('');

# Medikamentenvolumen auf gesamten Matrix
h_Matrix = 20000;
ri = 3200;
ro = 3400;
area_Matrix = 2*pi*h_Matrix*(ro+ri);
disp('Das kummulierte Medikamentenvolumen auf der gesamten Matrix beträgt ...[mu m^3].');
V_Matrix = V_REM*(area_Matrix/area_REM_mu2);
disp(V_Matrix)
disp('');

# Medikamentenmasse auf Matrix
rho_GM = 1.3*10^(-9);
m_Matrixoberflaeche = V_Matrix*rho_GM;
disp('Es befindet sich eine Medikamtenmasse von ...[mg] auf der Matrixoberfläche.');
disp(m_Matrixoberflaeche)
disp('');

# prozentualer Anteil der gesamten Medikamentenmasse auf der Matrixoberflaeche
m_d = 7.5;
q_GM_Oberflaeche = (m_Matrixoberflaeche/m_d)*100;
disp('Es befinden sich ...% des gesamten Medikaments auf der Matrixoberfläche.');
disp(q_GM_Oberflaeche)

# Gesamtanzahl an Medikamentenpartikeln auf der aeusseren Mantelflaeche
A_o = 2*pi*ro*h_Matrix;
n_o = A_o/area_REM_mu2 * anzahlPartikel;
disp('Auf der aeusseren Mantelflaeche befinden sich ... Partikel.');
disp(n_o)

# Gesamtanzahl an Medikamentenpartikeln auf der inneren Mantelflaeche
A_i = 2*pi*ri*h_Matrix;
n_i = A_i/area_REM_mu2 * anzahlPartikel;
disp('Auf der inneren Mantelflaeche befinden sich ... Partikel.');
disp(n_i)
