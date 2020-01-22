\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{IJFV}[2002/12/15 (F. Benkhaldoun and A. Rahmouni P13))]  
                 \LoadClass[titlepage,11pt]{article}
\ProcessOptions  
\RequirePackage{theorem}
\RequirePackage{amsmath}
\RequirePackage{psfrag}
\RequirePackage{fancyhdr}
\RequirePackage{amsfonts}
\RequirePackage{graphicx}
\theoremheaderfont{\scshape}
\theorembodyfont{\upshape}
\newcommand{\clearemptydoublepage}{%
        \newpage{\pagestyle{empty}\cleardoublepage}}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{definition}[theorem]{Definition\rm}
\newtheorem{remark}{\it Remark\/}
\newtheorem{example}{\it Example\/}
\renewcommand{\theequation}{\arabic{equation}}
\newcommand{\IJFVabstract}[1]{
\vspace{10mm}
\begin{center}
\begin{minipage}{.8\linewidth} 
{\bf Abstract}%\\[-3mm]

\vspace{2.5mm}
\hrule
\vspace{2mm}
{\small #1}%\\[-2mm]

\vspace{2mm}
\end{minipage}
\end{center}
}
\newcommand{\IJFVkeywords}[1]{
\begin{center}
\begin{minipage}{.8\linewidth}
{\small {\bf Key words} : #1}%\\[-2mm]

\vspace{2.5mm}
\hrule
\end{minipage}
\end{center}
\vspace{10mm}
}
\newcommand{\lag}{{\cal L}} 
\newcommand {\QED}{\hfill$\square$} 
\newcommand {\Cc}  {\mathbb{C}} 
\newcommand {\Qq}  {\mathbb{Q}} 
\newcommand {\Rr}  {\mathbb{R}} 
\newcommand {\Nn}  {\mathbb{N}}
\newcommand {\Pp}  {\mathbb{P}} 
\newcommand {\Zz}  {\mathbb{Z}} 
\newenvironment{proof}{{\bf Proof }}{\QED} 
\newcommand{\IJFVtitle}[2]{
\begin{center}
\thispagestyle{empty}
{\bf \Large #2}
\vspace{10mm}
\end{center}
}
\let\MyTitle\IJFVtitle
\renewcommand{\IJFVtitle}[1]{\newcommand{\IJFVotitle}[1]{#1}\MyTitle{#1}} 

%%%%%%% En tete etc %%%%%%%%%

\lhead[\it\leftmark]{}
\rhead[]{\small \IJFVotitle}
\rfoot[]{\thepage}
\cfoot[]{}
\lfoot[]{\small International Journal on Finite Volumes}  
\renewcommand {\footrulewidth} {0.4pt}   
\renewcommand {\plainheadrulewidth} {0.4pt}   
\renewcommand {\plainfootrulewidth} {0.4pt} 


\newcommand{\IJFVauthor}[1]{
\begin{center}
\thispagestyle{empty}
{\rm #1}\\
\vspace*{-3mm}
\end{center}
}
\newcommand{\IJFVinstitution}[1]{
\begin{center}
\thispagestyle{empty}
{\it \small #1}
\vspace*{-3mm}
\end{center}
}
\newcommand{\IJFVemail}[1]{
\begin{center}
\thispagestyle{empty}
{\small #1}\\
\vspace*{4mm}
\end{center}
}
\textheight 22cm
\textwidth 14cm
\pagestyle{fancyplain}
