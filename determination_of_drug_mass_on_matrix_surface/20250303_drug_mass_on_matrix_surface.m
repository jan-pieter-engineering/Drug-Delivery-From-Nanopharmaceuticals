clc
close all
clear all
pkg load image

% Read REM picture
REM=imread('REM_ohne_Legende.png');

% --------------------------------------------
% Convert REM picture into processible format
% --------------------------------------------
% 1. Convert into grayscale picture
REM_gray=rgb2gray(REM);
figure,
imshow(REM_gray);

% 2. Convert into black and white picture
threshold=0.78;
REM_bw=im2bw(REM_gray, threshold);
figure,
imshow(REM_bw);


% Count number of particles on REM picture [-]
anzahl = bwconncomp(REM_bw);
anzahlPartikel = anzahl.NumObjects;
disp('Die Anzahl der einzelnen Medikamentenpartikel auf der REM-Aufnahme betraegt:');
disp(anzahlPartikel)
disp('');
disp('');


% --------------------------------------------
% Specify REM area
% --------------------------------------------
disp('R E M  -  A U F N A H M E');

% 1. Calc total area of particles [pixels]
area = bwarea(REM_bw);
disp('Die Flaeche, die alle auf dem Ausschnitt zu sehenden Partikel zusammen einnehmen, betraegt ... [Pixel]');
disp(area);
display('');

% 2. Calc total area of REM picture [pixels]
area_REM_Pixel = 1072*730;
area_REM_mu2 = 221.58;
disp('Die REM-Aufnahme hat eine Auflösung von 1072x730 Pixeln und besteht damit aus ... [Pixel]');
disp(area_REM_Pixel)
display('');

% 3. Calc total area of particles [mu m^2]
disp('Die Fläche, die alle auf dem Ausschnitt zu sehenden Partikel zusammen einnehmen, betraegt ... [mu m^2]');
area_mu2 = (area_REM_mu2/area_REM_Pixel)*area;
disp(area_mu2)
disp('');
disp('');


% --------------------------------------------
% Specify particle geometry
% --------------------------------------------
disp('P A R T I K E L G E O M E T R I E');

% 1. Calc average particle diameter [mu m²]
area_mu2_pro_partikel = area_mu2/anzahlPartikel;
disp('Der durchschnittliche Partikelquerschnitt beträgt ... [mu m²].');
disp(area_mu2_pro_partikel)
disp('');

% 2. Calc portion of particles in REM picture [%]
area_prozent = area/area_REM_Pixel*100;
disp('Damit nehmen die Medikamentenpartikel ... [%] der Oberflaeche ein.');
disp(area_prozent)
disp('');

% 3. Calc average particle raddius [mu m]
disp('Der durchschnittliche Partikelradius beträgt ... [mu m].');
r = sqrt(area_mu2_pro_partikel/pi);
disp(r)
disp('');

% 4. Calc average particle surface [mu m^2]
disp('Die durchschnittliche Partikeloberfläche beträgt ... [mu m^2].');
a = 4*pi*r^2;
disp(a)
disp('');
disp('');


% --------------------------------------------
% Specify drug loading
% --------------------------------------------
disp('M E D I K A M E N T E N L A D U N G');

% 1. Calc accumulated drug volume on REM picture [mu m^3]
disp('Das kumulierte Medikamentenvolumen auf der REM-Aufnahme beträgt ... [mu m^3].');
V_REM = (4/3)*pi*r^3*anzahlPartikel;
disp(V_REM)
disp('');

% 2. Calc accumulated drug volume on whole matrix [mu m^3]
h_Matrix = 20000;
ri = 3200;
ro = 3400;
area_Matrix = 2*pi*h_Matrix*(ro+ri);
disp('Das kummulierte Medikamentenvolumen auf der gesamten Matrix beträgt ... [mu m^3].');
V_Matrix = V_REM*(area_Matrix/area_REM_mu2);
disp(V_Matrix)
disp('');

% 3. Calc drug mass on whole matrix [mg]
rho_GM = 1.3*10^(-9);
m_Matrixoberflaeche = V_Matrix*rho_GM;
disp('Es befindet sich eine Medikamtenmasse von ... [mg] auf der Matrixoberfläche.');
disp(m_Matrixoberflaeche)
disp('');

% 4. Calc drug portion on the matrix surface [%]
m_d = 7.5;
q_GM_Oberflaeche = (m_Matrixoberflaeche/m_d)*100;
disp('Es befinden sich ... [%] des gesamten Medikaments auf der Matrixoberfläche.');
disp(q_GM_Oberflaeche)
disp('');
