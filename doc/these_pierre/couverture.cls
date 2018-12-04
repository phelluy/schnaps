\NeedsTeXFormat{LaTeX2e}
 
\ProvidesClass{Couv}[23/04/2010 Classe personnelle, V1.2]
 
% classe de base
\LoadClass[a4paper,12pt,table]{book}
 
% extensions necessaires
\usepackage[utf8]{inputenc}
\RequirePackage[french]{babel}
\usepackage[utf8]{inputenc}
\usepackage[french]{babel}

%\RequirePackage{wallpaper}
\RequirePackage{anysize}
\RequirePackage[absolute]{textpos}
\RequirePackage{color}
\RequirePackage{tikz} 
\usetikzlibrary{calc}

% definition des textes
\def\nom{INSTITUT DE RECHERCHE MATH\'EMATIQUE AVANC\'EE}
\def\umr{UMR 7501}
\def\ville{Strasbourg}
\def\adresse{
	Universit\'e de Strasbourg et CNRS \\
 	7 Rue Ren\'e Descartes \\
 	67 084 STRASBOURG CEDEX
}
\def\adresseweb{www-irma.u-strasbg.fr}
\def\tel{T\'el. 03 68 85 01 29}
\def\fax{Fax 03 68 85 03 28}
\def\courriel{irma@math.unistra.fr}
\def\issn{ISSN 0755-3390}

% definition des parametres de la page de garde (en centimetres)
% definitions generales 
\def\margeint{2.2}	% marge de interieure
\def\margeext{1}	% mage de exterieure
\def\margehaut{1}	% marge superieure
\def\margebas{1}		% marge inferieure
% defintions specifiques pour la premiere de couverture
\def\hauteurlogo{3}		% hauteur de l'espace reserve aux logos
\def\largeurbande{5.2}	% largeur de la bande bleue foncee
% definitions specifiques pour la quatrieme de couverture
\def\hauteurcadre{12.5}	% hauteur du cadre avec le rappel des references du labo

\makeatletter
  
% reglages du positionnement absolu
\setlength{\TPHorizModule}{10mm} 
\setlength{\TPVertModule}{\TPHorizModule} 
\textblockorigin{0mm}{0mm} 

% definitions de couleurs
\definecolor{bleufonce}{rgb}{0.11,0.38,0.72}

% conversion en longueur pour les besoins latex
\newlength{\margeintlength}
\newlength{\margeextlength}
\newlength{\margehautlength}
\newlength{\margebaslength}
\newlength{\hauteurlogolength}
\newlength{\largeurbandelength}
\newlength{\hauteurcadrelength}

\setlength{\margeintlength}{\margeint cm}
\setlength{\margeextlength}{\margeext cm}
\setlength{\margehautlength}{\margehaut cm}
\setlength{\margebaslength}{\margebas cm}
\setlength{\hauteurlogolength}{\hauteurlogo cm}
\setlength{\largeurbandelength}{\largeurbande cm}
\setlength{\hauteurcadrelength}{\hauteurcadre cm}

% premiere de couverture
% calcul de la largeur de texte restante dans le cadre bleu clair
% on prevoit une marge de 2mm de chaque cote pour le texte situe dans le cadre bleu clair
\newlength{\largeurtextelength}
\setlength{\largeurtextelength}{20.6 cm}
\addtolength{\largeurtextelength}{-\margeintlength}
\addtolength{\largeurtextelength}{-\margeextlength}
\addtolength{\largeurtextelength}{-\largeurbandelength}

% quatrieme de couverture
% calcul de la largeur de texte disponible pour le resume
\newlength{\largeurresumelength}
\setlength{\largeurresumelength}{21 cm}
\addtolength{\largeurresumelength}{-\margeintlength}
\addtolength{\largeurresumelength}{-\margeextlength}

% definition des commandes specifiques : ury, resume, numero d'ordre et adresse tel
\newcommand{\jury}[1]{\def\@jury{#1}}
\newcommand{\resume}[1]{\def\@resume{#1}}
\newcommand{\issnordre}[1]{\def\@issnordre{#1}}
\newcommand{\telnum}[1]{\def\@telnum{#1}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{fancyhdr}
\RequirePackage{tabularx}
\RequirePackage{ifthen}
\RequirePackage[includefoot,nomarginpar,twoside,
    top=27mm,bottom=20mm,
    head=5mm,headsep=7mm,
    footskip=10mm,
    hmargin=20mm,bindingoffset=10mm]{geometry}

\newif\if@blankemptypage
\DeclareOption{noblankemptypage}{\@blankemptypagefalse}
\DeclareOption{blankemptypage}{\@blankemptypagetrue}

\@blankemptypagefalse
%\ExecuteOptions{}
\ProcessOptions

% names
\newcommand{\@titleapp}{Titre}
\newcommand{\@engtitleapp}{Title}
\newcommand{\@resumeapp}{R\'esum\'e}
\newcommand{\@abstractapp}{Abstract}
\newcommand{\@keywordsapp}{Mot-clefs}
\newcommand{\@engkeywordsapp}{Keywords}
\newcommand{\@juryapp}{Jury}
\newcommand{\@advisorapp}{Directeur de th\`ese}
\newcommand{\@coadvisorapp}{Co-directeur de th\`ese}
\newcommand{\@labapp}{Laboratoire}

% Sets the name of the laboratory
\newcommand{\@labo}{\textbf{!\texttt{labo} \`a d\'efinir!}}
\newcommand{\labo}[1]{\renewcommand{\@labo}{#1}}

% Sets the name of the school
\newcommand{\@school}{\textbf{!\texttt{school} \`a d\'efinir!}}
\newcommand{\school}[1]{\renewcommand{\@school}{#1}}

% Sets the name of the phd speciality
\newcommand{\@speciality}{\textbf{!\texttt{speciality} \`a d\'efinir!}}
\newcommand{\speciality}[1]{\renewcommand{\@speciality}{#1}}

% Sets the name of the university
\newcommand{\@universityabbrv}{\textbf{!\texttt{universit\'e} \`a d\'efinir!}}
\newcommand{\@university}{\textbf{!\texttt{universit'e} \`a d\'efinir!}}
\newcommand{\university}[2]{
    \renewcommand{\@universityabbrv}{#1}
    \renewcommand{\@university}{#2}
}

% Sets the ISBN number (if not set, prints lines for space to the 10 digits to
% be written
\newlength{\@ISBNcolwidth}
\setlength{\@ISBNcolwidth}{.25em}
\newcommand{\@ISBN}{
    \begin{tabular}{*{13}{|p{\@ISBNcolwidth}}|}
        &&&&&&&&&\\
        \hline
    \end{tabular}
}
\newcommand{\ISBN}[1]{\renewcommand{\@ISBN}{\texttt{#1}}}

% Sets the advisor name (and title, optional, defaults to M)
\newcommand{\@advisor}{\textbf{!\texttt{advisor} \`a d\'efinir!}}
\newcommand{\@advisortitle}{M}
\newcommand{\advisor}[2][M]{
    \renewcommand{\@advisortitle}{#1}
    \renewcommand{\@advisor}{#2}
}
	
% Sets the coadvisor name (and title, optional, defaults to M)
\newcommand{\@coadvisor}{\textbf{!\texttt{coadvisor} \`a d\'efinir!}}
\newcommand{\@coadvisortitle}{M}
\newcommand{\coadvisor}[2][M]{
    \renewcommand{\@coadvisortitle}{#1}
    \renewcommand{\@coadvisor}{#2}
}

% Defines a member of the jury
\newcommand{\jurymember}[3][M]{#1. &#2, &#3}
\newcommand{\juryadvisor}[1][Directeur de th\`ese]{
    \jurymember[\@advisortitle]{\@advisor}{#1}
}
\newcommand{\jurycoadvisor}[1][Co-directeur de th\`ese]{
    \jurymember[\@coadvisortitle]{\@coadvisor}{#1}
}

\newcommand{\@labaddr}{}
\newcommand{\labaddr}[1]{
    \renewcommand{\@labaddr}{#1}
}
\newcommand{\@engtitle}{}
\newcommand{\engtitle}[1]{
    \renewcommand{\@engtitle}{#1}
}


\newcommand{\@keywords}{}
\newcommand{\@engkeywords}{}
\newcommand{\keywords}[2]{
    \renewcommand{\@keywords}{#1}
    \renewcommand{\@engkeywords}{#2}
}

\newcommand{\@resumesize}{\small}
\newcommand{\resumesize}[1]{\renewcommand{\@resumesize}{#1}}
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% realisation de la premiere de couverture
% tronc commun : bandeau bleu, IRMA, logos, etc
\newcommand{\couvfond}{
\begin{tikzpicture}[remember picture,overlay] 
	% definition des coordonnees
	\draw (current page.north west) +(\margeint,-\margehaut) coordinate (A);	% coin superieur gauche du cadre bleu clair
	\draw (current page.north east) +(-\margeext -\largeurbande, -\margehaut) coordinate (B);	% coin superieur gauche du cadre bleu fonce
	\draw (current page.south east) +(-\margeext, \margebas + \hauteurlogo) coordinate (C);		% coin inferieur droit des deux cadres bleus
	\draw (current page.north east) +(-\margeext - 0.5*\largeurbande,-\margehaut) coordinate (I);	% point de depart pour les references de l'insitut
	\draw (current page.south east) +(-\margeext - 0.5*\largeurbande, \margebas+\hauteurlogo) coordinate (J);	% point de depart pour l'adresse web de l'institut.

	% dessin du cadre bleu clair
	\fill[color=bleufonce, opacity=0.1] (A) rectangle (C);
	% dessin du cadre bleu foncé
	\fill[color=bleufonce] (B) rectangle (C);	
	
	% mise en place des logos
	% logo de l'UDS
	\draw (current page.south east) +(-\margeext -\largeurbande ,\margebas+ 0.5*\hauteurlogo) node[anchor=center] (UDS) {\includegraphics[height=0.8\hauteurlogolength]{fig/cover/logoUDS.pdf}};
	% logo du CNRS
	\draw (current page.south west) +(\margeint ,\margebas+ 0.5*\hauteurlogo) node[anchor=west] (CNRS) {\includegraphics[height=0.8\hauteurlogolength]{fig/cover/logoCNRS.jpg}};
	% logo du CNRS (a mi-chemin entre le logo du CNRS et celui de l'UDS)
	\draw ($0.5*(CNRS.east)+0.5*(UDS.west)$) node[anchor=center]{\includegraphics[height=0.8\hauteurlogolength]{fig/cover/logoIRMA.png}};
	
	% mise en place des textes de la bande bleue foncee
	% nom et references de l'institut
	\draw (I) +(0,-0.3) node[anchor=north, text=white, text width=0.9\largeurbandelength, text badly centered]{
		\fontfamily{phv}
		\fontsize{16}{23}
		\selectfont
		\nom
		\\ \vspace*{6mm} 
		\umr 
		\\ \vspace*{6mm} 
		\ville
	};
	% adresse web
	\draw (J) +(0,+0.3) node[anchor=south, text=white, text width=0.9\largeurbandelength, text badly centered]{
		\fontfamily{phv}
		\fontsize{13}{17}
		\selectfont
		\adresseweb
	};
\end{tikzpicture} 
}

% Mise en place du texte principal pour les theses et les HDR
\newcommand{\couvtexte}{
	\begin{tikzpicture}[remember picture,overlay] 
		\draw (B) +(-0.2,0) node[anchor=north east, text=bleufonce, text width=\largeurtextelength]{
			\fontfamily{phv}%
			\fontsize{12}{18}
			\selectfont%
			\raggedleft%
			\vspace{5cm}
			{\bfseries \@publitype} \\
			\vspace*{5mm}
			{\large \@publichapeau} \\
			{\large \@specialite} \\
			\vspace*{8mm}
			{\large \bfseries \@author} \\
			\vspace{1cm}
			{\large \bfseries \@title} \\
			\vspace{4cm}
			{\large \@datecomplete}\\
			\vspace{1cm}
			\@jury\\
		};
	\end{tikzpicture}
}


% realisaton de la premiere de couverture
% cas d'une these
% specialite maths fonda

\newcommand{\makethesemath}{
	\def\@publitype{\huge Th\`ese}
	\def\@publichapeau{pr\'esent\'ee pour obtenir le grade de docteur de l'Universit\'e de Strasbourg}
	\def\@specialite{Sp\'ecialit\'e MATH\'EMATIQUES APPLIQU\'EES}
	\def\@datecomplete{Soutenue le \@date\\ devant la commission d'examen}

	\newpage
	\thispagestyle{empty}
	\couvfond %
	\couvtexte %
}

% realisaton de la premiere de couverture
% cas d'une these
% specialite maths app

\newcommand{\makethesemathapp}{
	\def\@publitype{\huge Th\`ese}
	\def\@publichapeau{pr\'esent\'ee pour obtenir le grade de docteur de l'Universit\'e de Strasbourg}
	\def\@specialite{Sp\'ecialit\'e MATH\'EMATIQUES APPLIQU\'EES}
	\def\@datecomplete{Soutenue le \@date\\ devant la commission d'examen}

	\newpage
	\thispagestyle{empty}
	\couvfond %
	\couvtexte %
}

% realisaton de la premiere de couverture
% cas d'une HDR
% specialite maths fonda

\newcommand{\makehdrmath}{
	\def\@publitype{\Large Habilitation \`a diriger des recherches}
	\def\@publichapeau{Universit\'e de Strasbourg}
	\def\@specialite{Sp\'ecialit\'e MATH\'EMATIQUES}
	\def\@datecomplete{Soutenue le \@date\\ devant la commission d'examen}

	\newpage
	\thispagestyle{empty}
	\couvfond %
	\couvtexte %
}

% realisaton de la premiere de couverture
% cas d'une HDR
% specialite maths app

\newcommand{\makehdrmathapp}{
	\def\@publitype{\Large Habilitation \`a diriger des recherches}
	\def\@publichapeau{Universit\'e de Strasbourg}
	\def\@specialite{Sp\'ecialit\'e MATH\'EMATIQUES APPLIQU\'EES}
	\def\@datecomplete{Soutenue le \@date\\ devant la commission d'examen}

	\newpage
	\thispagestyle{empty}
	\couvfond %
	\couvtexte %
}

% realisaton de la premiere de couverture
% cas d'une prepublication

\newcommand{\makeprepubli}{
	\newpage
	\thispagestyle{empty}
	\couvfond %
	\begin{tikzpicture}[remember picture,overlay] 
		\draw (B) +(-0.2,0) node[anchor=north east, text=bleufonce, text width=\largeurtextelength]{
			\fontfamily{phv}%
			\fontsize{12}{18}
			\selectfont%
			\raggedleft%
			\vspace{5cm}
			{\huge \bfseries Pr\'epublication} \\
			\vspace*{5mm}
			\vspace*{8mm}
			{\large \bfseries \@author} \\
			\vspace{1cm}
			{\large \@title} \\
			\vspace{4cm}
			{\large \@date}\\
		};
	\end{tikzpicture}
}


% realisation de la quatrieme de couverture
\newcommand{\makeresume}{
	\newpage
	\thispagestyle{empty}
	\begin{tikzpicture}[remember picture,overlay] 
		% coordonnees des points principaux
		\draw (current page.south west) +(\margeext,\margebas) coordinate (U);	% coin inferieur gauche du cadre
		\draw (current page.south east) +(-\margeint,\margebas+\hauteurcadre) coordinate (V);	% coin superieur gauche du cadre
		\draw (current page.south east) +(-\margeint,\margebas) coordinate (W);	% coin inferieur gauche du cadre
		\draw (current page.north west) +(\margeext,-\margehaut) coordinate (O);	% coin supereur gauche du resume
		
		% resume de la these
		\draw (O) node[anchor=north west, text width=\largeurresumelength, text justified] {
			\fontfamily{phv}%
			\fontsize{13}{18}
			\selectfont%
			\@resume\\
		};
		
		% mise en place du fond du cardre avec les formules mathematiques
		\begin{scope}
			\clip (U) rectangle (V);
			\draw (current page.south west) node[anchor=south west]  (fondA) {\includegraphics{fig/cover/formules.png}};
			\draw (fondA.north west) +(0,-0.2) node[anchor=south west]  (fondB) {\includegraphics{fig/cover/formules.png}};
			\draw (fondB.north west) +(0,-0.2) node[anchor=south west]  (fondC) {\includegraphics{fig/cover/formules.png}};
			\draw (fondC.north west) +(0,-0.2) node[anchor=south west]  (fondD) {\includegraphics{fig/cover/formules.png}};
		\end{scope}
		% dessin du cadre rectangulaire
		\draw[line width=1.5pt, color=bleufonce] (U) rectangle (V);
		
		% adresse du labo
		\draw (V) +(-0.2,-0.2) node[anchor=north east, text width=\largeurresumelength, color=bleufonce] {
			\fontfamily{phv}%
			\fontsize{13}{18}
			\selectfont%
			\raggedleft
			\nom\\
			\umr\\
			\adresse\\
			\vspace{5mm}
			\tel\\
			\fax\\
			\vspace{5mm}
			\adresseweb\\
			\courriel\\
		};
		
		% adresse tel et ordre dans la serie ISSN
		\draw (W) +(-0.2,0.2) node[anchor=south east, text width=\largeurresumelength, color=bleufonce] {
			\fontfamily{phv}%
			\fontsize{13}{18}
			\selectfont%
			\raggedleft
			IRMA \@issnordre \\
			\@telnum \\
		};
		
		% numero ISSN
		\draw (U) +(0.2,0.2) node[anchor=south west, text width=\largeurresumelength, color=bleufonce] (ISSN) {
			\fontfamily{phv}%
			\fontsize{13}{18}
			\selectfont%
			\noindent \issn
		};
		
		% logos
		% logo IRMA
		\draw (ISSN.north west) +(0,0.8) node[anchor=south west] (IRMA) {\includegraphics[width=3cm]{fig/cover/logoIRMA.pdf}};
		% logo UDS
		\draw (IRMA.north west) +(0,0.8) node[anchor=south west] (UDS) {\includegraphics[width=3cm]{fig/cover/logoUDS.pdf}};
		% logo CNRS
		\draw (UDS.north west) +(0,0.8) node[anchor=south west] (CNRS) {\includegraphics[width=3cm]{fig/cover/logoCNRS.pdf}};
		
		
	\end{tikzpicture}
}


% Back page
\newcommand{\@backsection}[4][r]{%
    \ifthenelse{\equal{#3}{}}{%
    }{%
        \ifthenelse{\equal{#2}{}}{}{\noindent\textbf{\textsc{#2}}%
            \ifthenelse{\equal{#1}{r}}{\\}{}}%
        \ifthenelse{\equal{#1}{r}}{\indent}{\noindent}{#4{#3}}\\%
        \ifthenelse{\equal{#1}{r}}{%
            \vspace*{\stretch{1}}%
            \noindent\rule{\hsize}{1pt}%
            \vspace*{\stretch{1}}%
        }{}%
     }%
}
\newcommand{\makeback}{
    \begin{titlepage}
        \thispagestyle{empty}
        \null\clearpage
        \thispagestyle{empty}
%        \vsize = 277mm
%        \voffset = -15mm
%        \topmargin = 0mm
%        \headheight = 0mm
%        \headsep = 0mm
%        \hsize = 160mm
%        \hoffset = -10mm
%        \vbox to \vsize {
            \noindent\rule{\hsize}{1pt}
            \vspace*{\stretch{1}}
            \@backsection{\@resumeapp}{\@resume}{\@resumesize}
            \@backsection{\@keywordsapp}{\@keywords}{\@resumesize}
            \@backsection{\@engtitleapp}{\@title}{\bf}
            \@backsection{\@abstractapp}{\@abstract}{\@resumesize}
            \@backsection{\@engkeywordsapp}{\@engkeywords}{\@resumesize}
            \@backsection[]{Adrr : }{\@labaddr}{\@resumesize}%
            \@backsection[]{ISBN : }{\@ISBN}{}
%            \vspace*{\stretch{1}}
%            \noindent\textbf{\textsc{ISBN}} : {\@resumesize{\@ISBN}}
%            \noindent\textbf{\textsc{ISBN}} : {\@resumesize{\@ISBN}}
    \end{titlepage}
}

% fancy pagestyle redefinition
    \fancyhead[LE,RO]{\small\leftmark}
    \fancyhead[LO,C,RE]{}
    \fancyfoot[LE,RO]{\rm\thepage}
    %\fancyfoot[LE]{\rm\thepage}
    \fancyfoot[LO,RE]{\small\rightmark}
    %\fancyfoot[LO,RE]{\hspace*{1cm}\small\rightmark}
    %\fancyfoot[LO]{\small\rightmark}
    \fancyfoot[C]{}
    %\renewcommand{\headrulewidth}{0.3pt}
    \renewcommand{\footrulewidth}{0.3pt}
% plain pagestyle redefinition
\fancypagestyle{plain}{
    \fancyhead[L,C,R]{}
    \fancyfoot[LE,RO]{\rm\thepage}
    \fancyfoot[LO,C,RE]{}
    \renewcommand{\headrulewidth}{0pt}
}
\pagestyle{fancy}

\renewcommand\thepart           {\@Roman\c@part}
\renewcommand\thechapter        {\@arabic\c@chapter}
\renewcommand\thesection        {\@arabic\c@section}
\renewcommand\thesubsection     {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection  {\thesubsection.\@arabic\c@subsubsection}
\renewcommand\theparagraph      {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph   {\theparagraph.\@arabic\c@subparagraph}

\newcommand{\openany}{\@openrightfalse}
\newcommand{\openright}{\@openrighttrue}

\newif{\if@frontmatter}
\renewcommand{\frontmatter}{
    \cleardoublepage
    \@mainmatterfalse
    \@frontmattertrue
    \@backmatterfalse
    \pagenumbering{roman}
    \fancyhead[LE,RO]{}
    \fancyhead[LO,C,RE]{}
    \fancyfoot[LE,RO]{}
    %\fancyfoot[LE]{\rm\thepage}
    \fancyfoot[LO,RE]{}
    %\fancyfoot[LO,RE]{\hspace*{1cm}\small\rightmark}
    %\fancyfoot[LO]{\small\rightmark}
    \fancyfoot[CE,CO]{\rm\thepage}

    \fancypagestyle{plain}{
      \fancyhead[L,C,R]{}
      %\fancyfoot[CO,CE]{\rm\thepage}
      \fancyfoot[C]{\rm\thepage}
      \renewcommand{\headrulewidth}{0pt}
}

    \renewcommand{\headrulewidth}{0.pt}
    \renewcommand{\footrulewidth}{0.pt}
}
\renewcommand{\mainmatter}{
    \cleardoublepage
    \@mainmattertrue
    \@frontmatterfalse
    \@backmatterfalse
    \pagenumbering{arabic}

    \fancyhead[LE,RO]{\small\leftmark}
    \fancyhead[LO,C,RE]{}
    \fancyfoot[LE,RO]{\rm\thepage}
    %\fancyfoot[LE]{\rm\thepage}
    \fancyfoot[LO,RE]{\small\rightmark}
    %\fancyfoot[LO,RE]{\hspace*{1cm}\small\rightmark}
    %\fancyfoot[LO]{\small\rightmark}
    \fancyfoot[C]{}
\fancypagestyle{plain}{
    \fancyhead[L,C,R]{}
    \fancyfoot[LE,RO]{\rm\thepage}
    \fancyfoot[LO,C,RE]{}
    \renewcommand{\headrulewidth}{0pt}
}

    \renewcommand{\headrulewidth}{0.3pt}
    \renewcommand{\footrulewidth}{0.3pt}
}
\newif{\if@backmatter}
\renewcommand{\backmatter}{
    \if@openright
        \cleardoublepage
    \else
        \clearpage
    \fi
    \@mainmatterfalse
    \@frontmatterfalse
    \@backmattertrue
}

\def\@chapter[#1]#2{
    \ifnum \c@secnumdepth >\m@ne
        \if@mainmatter
            \refstepcounter{chapter}
            \typeout{\@chapapp\space\thechapter.}
            \addcontentsline{toc}{chapter}
            {\protect\numberline{\thechapter}#1}
        \else
            \if@frontmatter
                \addcontentsline{toc}{section}{#1}
            \else
                \addcontentsline{toc}{chapter}{toto#1}
            \fi
        \fi
    \else
        \addcontentsline{toc}{chapter}{#1}
    \fi
    \chaptermark{#1}
    \addtocontents{lof}{\protect\addvspace{10\p@}}
    \addtocontents{lot}{\protect\addvspace{10\p@}}
    \if@twocolumn
        \@topnewpage[\@makechapterhead{#2}]
    \else
        \@makechapterhead{#2}
        \@afterheading
    \fi
}

\let\@oldschapter\@schapter
\def\@schapter#1{%
    \@oldschapter{#1}%
    \if@mainmatter
        %\addcontentsline{toc}{chapter}{{#1}}
    \else
        \if@frontmatter
            \addcontentsline{toc}{section}{#1}
        \else
           \if@backmatter
              \addcontentsline{toc}{part}{#1}
            \else
            \addcontentsline{toc}{chapter}{#1}
            \fi
        \fi
    \fi
    \@mkboth{\MakeUppercase{#1}}{\MakeUppercase{#1}}
}



\let\origdoublepage\cleardoublepage
\newcommand{\clearemptydoublepage}{%
\clearpage
{\pagestyle{empty}\origdoublepage}%
}
\let\cleardoublepage\clearemptydoublepage


\if@blankemptypage
    %Redefine cleardoublepage so that the pages inserted are really empty
    \renewcommand{\cleardoublepage}{
        \clearpage
        \if@twoside
            \ifodd
                \c@page
            \else
                \null
                \thispagestyle{empty} %set style empty
                \newpage
                \if@twocolumn\null\newpage\fi
            \fi
        \fi
    }
\fi



\makeatother




