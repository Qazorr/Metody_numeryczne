# Autor: Kacper Piątkowski #

## Pliki w katalogu:

#

* rozwiazania w pythonie do zadan widocznych nizej
* do zadan wymagajacych rysunki sporzadzone zostaly wykresy (pliki .ipynb)
* zdjecia przedstawiajace output programow

#

## 1 - rozwiązywanie układu równań ##

**(a) macierz trojdiagonalna:** \
\
$\begin{bmatrix}
4 & 1 & 0 & 0 & 0 & 0 & 0\\
1 & 4 & 1 & 0 & 0 & 0 & 0\\
0 & 1 & 4 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & 4 & 1 & 0 & 0\\
0 & 0 & 0 & 1 & 4 & 1 & 0\\
0 & 0 & 0 & 1 & 1 & 4 & 1\\
0 & 0 & 0 & 1 & 0 & 1 & 4
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3 \\
x_4 \\
x_5 \\
x_6 \\
x_7 
\end{bmatrix} =
\begin{bmatrix}
1 \\
2 \\
3 \\
4 \\
5 \\
6 \\
7
\end{bmatrix}$ \
\
Układ równan rozwiazany za pomocą [Algorytmu Thomasa](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm) \
Dzieki temu algorytmowi nasza zlozonosc wynosi $\mathcal{O}(n)$

**(b) taka sama macierz tylko z 1 w rogach** \
\
$\begin{bmatrix}
4 & 1 & 0 & 0 & 0 & 0 & 1\\
1 & 4 & 1 & 0 & 0 & 0 & 0\\
0 & 1 & 4 & 1 & 0 & 0 & 0\\
0 & 0 & 1 & 4 & 1 & 0 & 0\\
0 & 0 & 0 & 1 & 4 & 1 & 0\\
0 & 0 & 0 & 1 & 1 & 4 & 1\\
1 & 0 & 0 & 1 & 0 & 1 & 4
\end{bmatrix}
\begin{bmatrix}
x_1 \\
x_2 \\
x_3 \\
x_4 \\
x_5 \\
x_6 \\
x_7 
\end{bmatrix} =
\begin{bmatrix}
1 \\
2 \\
3 \\
4 \\
5 \\
6 \\
7
\end{bmatrix}$ \
\
Uzyskujemy macierz trojdiagonalna poprzez odjecie macierzy majacej 1 w kazdym z rogów \
Korzystamy z algorytmu [Shermana-Morissona](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad03.pdf) (str 55)

## 2 - rozwiazywanie ukladu rownan 128x128 ##
\
$\begin{pmatrix}
4 & 1 & 0 & 0 & 1 & \cdots & \cdots\\
1 & 4 & 1 & 0 & 0 & \cdots & \cdots\\
0 & 1 & 4 & 0 & 0 & \cdots & \cdots\\
0 & 0 & 1 & 4 & 1 & \cdots & \cdots\\
1 & 0 & 0 & 1 & 4 & \cdots & \cdots\\
\vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \ddots\\
\cdots & \cdots & 1 & 0 & 0 & 1 & 4
\end{pmatrix}$

* W zadaniu zostaly zaimplementowane algorytmy z wykladu 4
**(a)** 
[metoda Gaussa-Seidela](http://th-www.if.uj.edu.pl/zfs/gora/metnum18/wyklad04.pdf) (strona 4) - uwzgledniona zostala budowa macierzy (mnozenia przez 0 zostaly pominiete) \
**(b)** 
[metoda gradientow sprzężonych](http://th-www.if.uj.edu.pl/zfs/gora/metnum18/wyklad04.pdf) (strona 15)

## 3 - znajdowanie 2 najwiekszych wartosci wlasnych i ich wektorow ##

\
$\large \begin{bmatrix}
\dfrac{19}{12} & \dfrac{13}{12} & \dfrac{5}{6} & \dfrac{5}{6} & \dfrac{13}{12} & \dfrac{-17}{12} \\[10pt]
\dfrac{13}{12} & \dfrac{13}{12} & \dfrac{5}{6} & \dfrac{5}{6} & \dfrac{-11}{12} & \dfrac{13}{12} \\[10pt]
\dfrac{5}{6} & \dfrac{5}{6} & \dfrac{5}{6} & \dfrac{-1}{6} & \dfrac{5}{6} & \dfrac{5}{6} \\[10pt]
\dfrac{5}{6} & \dfrac{5}{6} & \dfrac{-1}{6} & \dfrac{5}{6} & \dfrac{5}{6} & \dfrac{5}{6} \\[10pt]
\dfrac{13}{11} & \dfrac{-11}{12} & \dfrac{5}{6} & \dfrac{5}{6} & \dfrac{13}{12} & \dfrac{13}{12} \\[10pt]
\dfrac{-17}{12} & \dfrac{13}{12} & \dfrac{5}{6} & \dfrac{5}{6} & \dfrac{13}{12} & \dfrac{19}{12}
\end{bmatrix}$

* W zadaniu nalezalo skorzystac z [metody potegowej](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad05.pdf) (strona 13), 
dodatkowo skorzystalem z [twierdzenia](https://math.stackexchange.com/questions/1114777/approximate-the-second-largest-eigenvalue-and-corresponding-eigenvector-given), 
dzieki ktoremu obliczylem druga wartosc wlasna

## 4 - sprowadzenie do postaci trojdiagonalnej i obliczenie wartosci wlasnych ##
* Posluzylem sie [transformacja householdera](https://johnfoster.pge.utexas.edu/numerical-methods-book/LinearAlgebra_EigenProblem2.html) 
dzieki ktorej sprowadzilem macierz do postaci trojdiagonalnej \
W nastepnej czesci skorzystalem z [Obrotów Givensa](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad03.pdf) (strona 45) \
oraz [algorytmu QR dla macierzy symetrycznych trojdiagonalnych](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad05.pdf) (strona 26) 

## 6  - znajdowanie wektora wlasnego majac wartosc wlasna ##
* Wartosc wlasna $\lambda \approx 0.38197$
* W zadaniu mamy podana macierz oraz wartosc wlasna dla ktorej mamy wyznaczyc wektor wlasny \
Skorzystalem z [wykladu 5](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad05.pdf) (strona 43/44) na ktorym zostal opisany algorytm 
[Inverse iteration](https://en.wikipedia.org/wiki/Inverse_iteration) dzieki ktoremu obliczylem wektor wlasny

## 7 - znajdowanie wielomianu interpolacyjnego opartego na tabeli ##


| x     | 0.062500 | 0.187500 |  0.312500 | 0.437500 | 0.562500 | 0.687500 | 0.812500 | 0.937500
| :---: |   :-:    |    :-:   |    :-:    | :-: | :-: | :-: | :-: | :-: |
| f(x)  | 0.687959 | 0.073443 | −0.517558 | −1.077264 | −1.600455 | −2.080815 | −2.507266 | −2.860307


* W zadaniu skonstruowalem [macierz Vandermonde'a](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf) 
(strona 14) z wezlow z tabeli po czym rozwiazalem uklad rownan, \
w ktorym wartosciami szukanymi byly wspolczynniki wielomianu

## 8 - konstruowanie wielomianu interpolacyjnego Lagrange'a ##

* W zadaniu mamy podane wezly oraz funkcje, dzieki ktorej obliczamy wartosci w naszych wezlach \
$\displaystyle f(x) = \frac{1}{1 + 5x^2}$ \
Nastepnie korzystam ze [wzoru interpolacyjnego Lagrange'a](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf) (strona 16) dzieki ktoremu otrzymuje wspolczynniki wielomianu

## 9 - konstruowanie splajnu kubicznego ##
* Tworzymy splajn kubiczny opierajac sie na wezlach z zadania `8.py` \
Korzystajac z [wykladu 6](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf) (strona 34) wyliczam **A, B, C, D** \
a nastepnie znajduje $\zeta$ z [rownania](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad06.pdf) (strona 38), poniewaz wiem ze nasze wezly sa rownoodalone od siebie ($\displaystyle \big|x_{i+1} - x_{i}\big| = \frac{2}{8}$) \
Na koniec wyliczam wielomiany dla kazdego odcinka miedzy kolejnymi punktami korzystajac ze wzoru ze strony 34.

## 11 - obliczanie calki metoda trapezow i metoda Romberga ##
* W zadaniu mamy zadana całke:
* $\displaystyle\int_{0}^{\infin} \sin(\pi\frac{1 + \sqrt{x}}{1 + x^2})e^{-x}$
* Zaczynam od wyliczenia **A** korzystajac z podpowiedzi w zadaniu \
Nastepnie korzystam ze [zlozonego wzoru trapezow](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad07.pdf) (strona 20) dzieki ktoremu wyliczam calke zageszczajac ilosc trapezow \
W drugiej czesci korzystam z [Metody Romberga](http://th-www.if.uj.edu.pl/zfs/gora/metnum21/wyklad07.pdf) (strony 28-30) sprawdzajac przekatna macierzy Romberga

## 12 - obliczenie lim x->∞ z calki od −∞ do x ## 
* W zadaniu mamy zadana funkcje:
$\newcommand{\fx}{cos(\frac{1+t}{t^2+0.04})e^{-t^2}}$
$\displaystyle f(t) = {cos(\frac{1+t}{t^2+0.04})e^{-t^2}}$ \
Nalezy policzyc: $\displaystyle\lim_{x \to \infin} F(x) = \int_{-\infty}^{x} {cos(\frac{1+t}{t^2+0.04})e^{-t^2}}$
1) Narysowalem funkcje i zauwazylem ze w $\displaystyle \lim_{x \to \pm\infin} f(x) \to 0$
2) Wyliczylem punkt startowy i koncowy dla calki szukajac punktu dla ktorego f(x) <= $10^{-8}$
3) Rozwiazuje calke metoda trapezow tak jak w zadaniu 11

* **[Side note]** \
Do tego samego mozemy dojsc analitycznie, poniewaz mnozymy przez $e^{-x^2}$ \
Wiemy ze $x^2$ jest zawsze dodatnie, wiec potega bedzie zawsze ujemna, \
co znaczy ze kolejne potegi beda dawaly nam mniejsze liczby \
Dlatego $-\infin$ jak i w $\infin$ nasza funkcja zbliza sie do zera, poniewaz $\displaystyle \lim_{x \to \pm\infin} e^{-x^2} \to 0$

## 13 - metoda Lagurre'a do rozwiazywania równań ##
* W zadaniu mamy 3 rownania do rozwiazania (1 i 2 posiadaja rzeczywiste wspolczynniki a 3 zespolone) \
1) $243z^7  − 486z^6 + 783z^5 − 990z^4 + 558z^3 − 28z^2 − 72z + 16 = 0$
2) $z^{10} + z^9 + 3z^8 + 2z^7 − z^6 − 3z^5 − 11z^4 − 8z^3 − 12z^2 − 4z − 4 = 0$
3) $z^4 + iz^3 − z^2 − iz + 1 = 0$

* Korzystam w zadaniu z [Metody Lagurre'a](http://th-www.if.uj.edu.pl/zfs/gora/metnum15/wyklad10.pdf) (strona 14) dzieki ktorej znajduje pojedyncze miejsca zerowe wielomianu \
Nastepnie dokonuje deflacje wielomianu (wyciagam znane juz miejsce zerowe z wielomianu) \
Po czym ponownie licze miejsce zerowe z Metody Lagurre'a \
Po wykonaniu wystarczajacej liczby takich iteracji otrzymuje wszystkie numeryczne miejsca zerowe wielomianu

## 15 - konstruowanie naturalnego splajnu kubiczny ##
* Wczytuje wezly z pliku, najpierw sprawdzilem czy punkty sa od siebie rownoodalone i ostatecznie \
korzystajac z rozwiazania zadania `9.py` tworze naturalny splajn kubiczny